version 1.0

##
## InspectGvcfs.wdl
##
## For each GVCF in the input list, reports:
##   1. variant_count   — number of variant records (reference blocks excluded)
##   2. sequencing_type — "WGS" or "WES", inferred from total genomic bases covered
##   3. is_reblocked    — true/false; detected via header metadata and block-size heuristic
##
## Per-sample outputs (one file per GVCF):
##   {sample_id}_inspect.txt — human-readable report
##   result.tsv              — machine-readable one-row TSV
##
## Workflow-level outputs:
##   per_sample_reports — Array[File] of *_inspect.txt, one per GVCF
##   per_sample_tsvs    — Array[File] of result.tsv rows, one per GVCF
##   summary_tsv        — all rows combined into one table
##
## Dependencies: bcftools (available in the default GATK docker image)
##

workflow InspectGvcfs {

  input {
    Array[File] gvcfs           # one .g.vcf.gz (or .g.vcf) per sample
    Array[File] gvcf_indices    # corresponding .tbi index files

    # Bases-covered threshold (Mb) above which a sample is classified as WGS.
    # WGS typically covers ~3,000 Mb; WES ~50–200 Mb.  Default: 1500 Mb.
    Int wgs_coverage_threshold_mb = 1500

    # Runtime
    String docker    = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    Int    cpu       = 1
    Int    memory_gb = 4
    Int    disk_gb   = 50   # increase for large WGS GVCFs (>5 GB each)
  }

  # ── Scatter: inspect each GVCF independently ────────────────────────────────
  Array[Pair[File, File]] gvcf_pairs = zip(gvcfs, gvcf_indices)

  scatter (pair in gvcf_pairs) {
    call InspectOneGvcf {
      input:
        gvcf                      = pair.left,
        gvcf_index                = pair.right,
        wgs_coverage_threshold_mb = wgs_coverage_threshold_mb,
        docker                    = docker,
        cpu                       = cpu,
        memory_gb                 = memory_gb,
        disk_gb                   = disk_gb
    }
  }

  # ── Gather: merge all per-sample TSV rows into one summary file ─────────────
  call GatherResults {
    input:
      per_sample_tsvs = InspectOneGvcf.result_tsv,
      docker          = docker,
      cpu             = 1,
      memory_gb       = 2,
      disk_gb         = 10
  }

  # ── Outputs ─────────────────────────────────────────────────────────────────
  output {
    # ── Per-sample files (one file per input GVCF, named after the sample) ───
    # Human-readable report: {sample_id}_inspect.txt
    Array[File]    per_sample_reports  = InspectOneGvcf.per_sample_report
    # Machine-readable one-row TSV: sample_id, variant_count, sequencing_type, is_reblocked
    Array[File]    per_sample_tsvs     = InspectOneGvcf.result_tsv

    # ── Per-sample scalar arrays (one element per input GVCF, in input order) ─
    Array[String]  sample_ids       = InspectOneGvcf.sample_id
    Array[Int]     variant_counts   = InspectOneGvcf.variant_count
    Array[String]  sequencing_types = InspectOneGvcf.sequencing_type  # "WGS" or "WES"
    Array[Boolean] is_reblocked     = InspectOneGvcf.is_reblocked

    # ── Aggregated summary (all samples in one file) ──────────────────────────
    # Columns: sample_id, variant_count, sequencing_type, is_reblocked
    File summary_tsv = GatherResults.summary_tsv
  }
}


# ────────────────────────────────────────────────────────────────────────────
task InspectOneGvcf {
# ────────────────────────────────────────────────────────────────────────────
  input {
    File   gvcf
    File   gvcf_index
    Int    wgs_coverage_threshold_mb
    String docker
    Int    cpu
    Int    memory_gb
    Int    disk_gb
  }

  command <<<
    set -euo pipefail
    GVCF="~{gvcf}"

    # ── 1. Sample ID from VCF header ─────────────────────────────────────────
    SAMPLE_ID=$(bcftools query -l "$GVCF" | head -1)
    echo "$SAMPLE_ID" > sample_id.txt

    # ── 2. Variant count ─────────────────────────────────────────────────────
    # Distinguish reference blocks from variant records by the ALT column (col 5):
    #   Reference block : ALT = <NON_REF>  (sole symbolic allele, starts with '<')
    #   Variant record  : ALT = A,<NON_REF> or A or AT,... (starts with a real base)
    #
    # NOTE: do NOT use `bcftools view -e 'ALT="<NON_REF>"'` here.
    # bcftools filter expressions treat '<' and '>' as comparison operators,
    # so the expression silently misfires — either keeping all records or none.
    # The awk column-5 check below is the reliable alternative.
    VAR_COUNT=$(bcftools view -H "$GVCF" \
      | awk 'BEGIN{c=0} $5 !~ /^</{c++} END{print c}')
    echo "$VAR_COUNT" > variant_count.txt

    # ── 3. WGS vs WES ────────────────────────────────────────────────────────
    # Sum total genomic bases covered:
    #   - Reference blocks carry INFO/END  → span = END − POS + 1
    #   - Variant records                  → span = length(REF)
    # An early-exit is used so large WGS files are not read to completion.
    THRESHOLD_BASES=$(( ~{wgs_coverage_threshold_mb} * 1000000 ))

    SEQ_TYPE=$(bcftools query -f '%POS\t%INFO/END\t%REF\n' "$GVCF" | \
      awk -v thr="$THRESHOLD_BASES" '
        BEGIN { seqtype = "WES"; bases = 0 }
        {
          if ($2 != ".") {
            bases += ($2 - $1 + 1)
          } else {
            bases += length($3)
          }
          if (bases > thr) { seqtype = "WGS"; exit }
        }
        END { print seqtype }
      ')
    echo "$SEQ_TYPE" > sequencing_type.txt

    # ── 4. Reblocked status ──────────────────────────────────────────────────
    # Primary check  : "ReblockGVCF" in any ##GATKCommandLine header line.
    # Secondary check: >10 % of reference blocks span >500 bp (reblocking merges
    #                  many small GQ bands into large blocks).
    IS_REBLOCKED="false"

    if bcftools view -h "$GVCF" | grep -qi "ReblockGVCF"; then
      IS_REBLOCKED="true"
    else
      IS_REBLOCKED=$(bcftools query -f '%POS\t%INFO/END\n' "$GVCF" | \
        awk '
          $2 != "." {
            total++
            if (($2 - $1 + 1) > 500) large++
          }
          END {
            if (total > 0 && large / total > 0.10) print "true"
            else                                    print "false"
          }
        ')
    fi
    echo "$IS_REBLOCKED" > is_reblocked.txt

    # ── Write one-row TSV for GatherResults ──────────────────────────────────
    printf '%s\t%s\t%s\t%s\n' \
      "$SAMPLE_ID" "$VAR_COUNT" "$SEQ_TYPE" "$IS_REBLOCKED" \
      > result.tsv

    # ── Write named human-readable report for this sample ────────────────────
    GVCF_BASENAME=$(basename "~{gvcf}")
    cat > "${SAMPLE_ID}_inspect.txt" << REPORT
=== GVCF Inspection Report ===
Sample ID       : ${SAMPLE_ID}
GVCF file       : ${GVCF_BASENAME}
Variant count   : ${VAR_COUNT}
Sequencing type : ${SEQ_TYPE}
Reblocked       : ${IS_REBLOCKED}
REPORT
  >>>

  runtime {
    docker: docker
    cpu:    cpu
    memory: "~{memory_gb} GiB"
    disks:  "local-disk ~{disk_gb} HDD"
  }

  output {
    String  sample_id         = read_string("sample_id.txt")
    Int     variant_count     = read_int("variant_count.txt")
    String  sequencing_type   = read_string("sequencing_type.txt")
    Boolean is_reblocked      = read_boolean("is_reblocked.txt")
    File    result_tsv        = "result.tsv"
    File    per_sample_report = glob("*_inspect.txt")[0]
  }
}


# ────────────────────────────────────────────────────────────────────────────
task GatherResults {
# ────────────────────────────────────────────────────────────────────────────
  input {
    Array[File] per_sample_tsvs
    String docker
    Int    cpu
    Int    memory_gb
    Int    disk_gb
  }

  command <<<
    set -euo pipefail
    printf 'sample_id\tvariant_count\tsequencing_type\tis_reblocked\n' > summary.tsv
    cat ~{sep=" " per_sample_tsvs} >> summary.tsv
  >>>

  runtime {
    docker: docker
    cpu:    cpu
    memory: "~{memory_gb} GiB"
    disks:  "local-disk ~{disk_gb} HDD"
  }

  output {
    File summary_tsv = "summary.tsv"
  }
}
