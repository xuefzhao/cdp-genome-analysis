version 1.0

## TargetDepthSummary
##
## Summarize WES sequencing depth inside vs. outside the capture targets.
##
## Input: a list of mosdepth per-base coverage BED files for ONE sample, each
## file representing one chromosome (e.g. <sample>.chr1.coverage.per-base.bed.gz),
## plus the WES target intervals BED. mosdepth per-base format is run-length
## encoded: columns  chrom  start  end  depth  (0-based half-open, integer depth).
##
## Output: mean / median / sd of per-base depth, computed over all bases
## (1) WITHIN the targets and (2) OUTSIDE the targets, aggregated across all the
## input chromosomes for the sample.
##
## Method: depth is base-weighted via a depth histogram (depth -> #bases), which
## is exact and memory-light even over whole-genome per-base files. Targets are
## sorted + merged so overlapping target intervals never double-count a base.
## "within"  = bedtools intersect (per-base clipped to targets)
## "outside" = bedtools subtract  (per-base minus targets)

workflow TargetDepthSummary {
  input {
    Array[File] per_base_beds          # one mosdepth per-base bed(.gz) per chromosome, single sample
    File target_bed                    # WES target intervals (BED or BED.gz)
    String sample_id = "sample"
    Int depth_col = 4                  # 1-based column holding depth in the per-base bed
    String bedtools_docker = "quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"
    String python_docker = "python:3.11-slim"
  }

  scatter (bed in per_base_beds) {
    call ChromHistograms {
      input:
        per_base_bed = bed,
        target_bed   = target_bed,
        depth_col    = depth_col,
        docker       = bedtools_docker
    }
  }

  call Summarize {
    input:
      within_hists  = ChromHistograms.within_hist,
      outside_hists = ChromHistograms.outside_hist,
      sample_id     = sample_id,
      docker        = python_docker
  }

  output {
    File summary_tsv  = Summarize.summary_tsv     # tidy table: region, n_bases, mean, median, sd
    Map[String, Float] summary = Summarize.summary

    # convenience scalar outputs
    Float within_mean    = Summarize.summary["within_mean"]
    Float within_median  = Summarize.summary["within_median"]
    Float within_sd      = Summarize.summary["within_sd"]
    Float outside_mean   = Summarize.summary["outside_mean"]
    Float outside_median = Summarize.summary["outside_median"]
    Float outside_sd     = Summarize.summary["outside_sd"]
  }
}

## Per-chromosome: split per-base depth into within/outside target and emit a
## depth histogram (depth<TAB>n_bases) for each.
task ChromHistograms {
  input {
    File per_base_bed
    File target_bed
    Int depth_col = 4
    String docker
    Int cpu = 1
    Int mem_gb = 4
  }
  Int disk_gb = ceil(size(per_base_bed, "GB") + size(target_bed, "GB")) * 5 + 20

  command <<<
    set -euo pipefail

    decompress() { case "$1" in *.gz) gzip -dc "$1";; *) cat "$1";; esac; }

    # Per-base coverage for this chromosome (sorted).
    decompress "~{per_base_bed}" | sort -k1,1 -k2,2n > perbase.bed

    # Targets: keep first 3 BED columns, sort, and MERGE so overlapping/adjacent
    # target intervals can't double-count bases in the intersect step.
    decompress "~{target_bed}" | awk -v OFS='\t' '$0!~/^(#|track|browser)/{print $1,$2,$3}' \
      | sort -k1,1 -k2,2n \
      | bedtools merge -i - > targets.merged.bed

    # WITHIN targets: intersection clips per-base intervals to target bounds; depth
    # retained in column ~{depth_col}. Weight depth by interval length (#bases).
    bedtools intersect -a perbase.bed -b targets.merged.bed \
      | awk -v OFS='\t' -v c=~{depth_col} '{h[$c]+=($3-$2)} END{for(d in h) print d,h[d]}' \
      | sort -k1,1n > within.hist

    # OUTSIDE targets: per-base intervals with the target portions removed.
    bedtools subtract -a perbase.bed -b targets.merged.bed \
      | awk -v OFS='\t' -v c=~{depth_col} '{h[$c]+=($3-$2)} END{for(d in h) print d,h[d]}' \
      | sort -k1,1n > outside.hist

    # Ensure non-empty files even if a region had no bases.
    touch within.hist outside.hist
  >>>

  output {
    File within_hist  = "within.hist"
    File outside_hist = "outside.hist"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: mem_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 2
  }
}

## Merge per-chromosome histograms and compute mean / median / sd of depth for
## the within-target and outside-target base populations.
task Summarize {
  input {
    Array[File] within_hists
    Array[File] outside_hists
    String sample_id = "sample"
    String docker
  }

  command <<<
    set -euo pipefail
    python3 <<'PY'
import json, math

within_list  = "~{write_lines(within_hists)}"
outside_list = "~{write_lines(outside_hists)}"
sample_id    = "~{sample_id}"

def load_hist(list_file):
    """Sum depth->n_bases across all per-chromosome histogram files."""
    h = {}
    with open(list_file) as lf:
        paths = [ln.strip() for ln in lf if ln.strip()]
    for p in paths:
        with open(p) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                depth_s, n_s = line.split("\t")
                d = int(depth_s)
                h[d] = h.get(d, 0) + int(n_s)
    return h

def stats(h):
    n = sum(h.values())
    if n == 0:
        return {"n_bases": 0, "mean": 0.0, "median": 0.0, "sd": 0.0}
    total = sum(d * c for d, c in h.items())
    mean = total / n
    var = sum(c * (d - mean) ** 2 for d, c in h.items()) / n   # population variance
    sd = math.sqrt(var)
    items = sorted(h.items())                                   # by depth
    def value_at(pos):              # pos is 1-based rank
        run = 0
        for d, c in items:
            run += c
            if run >= pos:
                return d
        return items[-1][0]
    if n % 2 == 1:
        median = float(value_at((n + 1) // 2))
    else:
        median = (value_at(n // 2) + value_at(n // 2 + 1)) / 2.0
    return {"n_bases": n, "mean": mean, "median": float(median), "sd": sd}

w = stats(load_hist(within_list))
o = stats(load_hist(outside_list))

summary = {
    "within_n_bases":  float(w["n_bases"]),
    "within_mean":     w["mean"],
    "within_median":   w["median"],
    "within_sd":       w["sd"],
    "outside_n_bases": float(o["n_bases"]),
    "outside_mean":    o["mean"],
    "outside_median":  o["median"],
    "outside_sd":      o["sd"],
}
with open("summary.json", "w") as f:
    json.dump(summary, f)

with open("summary.tsv", "w") as f:
    f.write("sample_id\tregion\tn_bases\tmean\tmedian\tsd\n")
    for region, s in [("within_targets", w), ("outside_targets", o)]:
        f.write(f"{sample_id}\t{region}\t{s['n_bases']}\t{s['mean']:.4f}\t{s['median']:.4f}\t{s['sd']:.4f}\n")

print(open("summary.tsv").read())
PY
  >>>

  output {
    File summary_tsv = "summary.tsv"
    Map[String, Float] summary = read_json("summary.json")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 20 HDD"
    preemptible: 2
  }
}
