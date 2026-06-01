version 1.0

##
## InspectGvcf.wdl
##
## Inspects a single GVCF and reports:
##   1. variant_count   — number of variant records (reference blocks excluded)
##   2. sequencing_type — "WGS" or "WES", inferred from total genomic bases covered
##   3. is_reblocked    — true/false; detected via header metadata and block-size heuristic
##
## Outputs:
##   sample_id, variant_count, sequencing_type, is_reblocked  (scalar values)
##   report_txt  — human-readable {sample_id}_inspect.txt
##   result_tsv  — machine-readable one-row TSV
##
## Dependencies: bcftools (available in the default GATK docker image)
##

workflow InspectGvcf {

  input {
    File gvcf        # single .g.vcf.gz (or .g.vcf)
    File gvcf_index  # corresponding .tbi index

    # Bases-covered threshold (Mb) above which a sample is classified as WGS.
    # WGS typically covers ~3,000 Mb; WES ~50–200 Mb.  Default: 1500 Mb.
    Int wgs_coverage_threshold_mb = 1500

    # Runtime
    String docker    = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    Int    cpu       = 1
    Int    memory_gb = 4
    Int    disk_gb   = 50   # increase for large WGS GVCFs (>5 GB each)
  }

  call InspectOneGvcf {
    input:
      gvcf                      = gvcf,
      gvcf_index                = gvcf_index,
      wgs_coverage_threshold_mb = wgs_coverage_threshold_mb,
      docker                    = docker,
      cpu                       = cpu,
      memory_gb                 = memory_gb,
      disk_gb                   = disk_gb
  }

  output {
    String  sample_id       = InspectOneGvcf.sample_id
    Int     variant_count   = InspectOneGvcf.variant_count
    String  sequencing_type = InspectOneGvcf.sequencing_type  # "WGS" or "WES"
    Boolean is_reblocked    = InspectOneGvcf.is_reblocked
    File    report_txt      = InspectOneGvcf.per_sample_report
    File    result_tsv      = InspectOneGvcf.result_tsv
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

    # set +o pipefail inside the subshell so that when awk exits early (WGS case),
    # bcftools gets SIGPIPE but that non-zero exit does NOT abort the outer script.
    SEQ_TYPE=$(set +o pipefail; bcftools query -f '%POS\t%INFO/END\t%REF\n' "$GVCF" | \
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

    # Read header into a variable first so grep has no upstream pipe to bcftools;
    # avoids SIGPIPE (bcftools exit 141) being mis-reported as a false condition
    # by pipefail.
    VCF_HEADER=$(bcftools view -h "$GVCF")
    if echo "$VCF_HEADER" | grep -qi "ReblockGVCF"; then
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

