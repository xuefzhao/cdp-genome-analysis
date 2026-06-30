version 1.0

## IntegrateMosdepthSummary
##
## Combine per-chromosome mosdepth `*.summary.txt` files for ONE sample into a
## single genome-wide summary, and recompute the genome mean coverage.
##
## Each mosdepth summary has columns:  chrom  length  bases  mean  min  max
## and contains the chromosome row plus a redundant `total` row. This workflow
## keeps the per-chromosome rows (dropping the redundant `total` rows), sums
## `length` and `bases` across chromosomes, and recomputes:
##     genome mean = sum(bases) / sum(length)        (mosdepth's own definition)
##     genome min  = min over chromosomes
##     genome max  = max over chromosomes
##
## Outputs:
##   1. genome_summary : one combined summary file (per-chrom rows + a recomputed `total` row)
##   2. mean           : the recalculated genome mean coverage, as a Float

workflow IntegrateMosdepthSummary {
  input {
    Array[File] per_chrom_summaries     # mosdepth .summary.txt, one per chromosome, single sample
    String sample_id
    String docker = "python:3.11-slim"
    Int cpu = 1
    Float mem_gb = 2.0
    Int disk_gb = 20
    Int preemptible_tries = 2
    Int max_retries = 1
  }

  call IntegrateSummary {
    input:
      per_chrom_summaries = per_chrom_summaries,
      sample_id = sample_id,
      docker = docker,
      cpu = cpu,
      mem_gb = mem_gb,
      disk_gb = disk_gb,
      preemptible_tries = preemptible_tries,
      max_retries = max_retries
  }

  output {
    File genome_summary = IntegrateSummary.genome_summary
    Float mean = IntegrateSummary.mean
  }
}

task IntegrateSummary {
  input {
    Array[File] per_chrom_summaries
    String sample_id
    String docker
    Int cpu
    Float mem_gb
    Int disk_gb
    Int preemptible_tries
    Int max_retries
  }

  command <<<
    set -euo pipefail
    python3 <<'PY'
list_file = "~{write_lines(per_chrom_summaries)}"
sample_id = "~{sample_id}"

paths = [ln.strip() for ln in open(list_file) if ln.strip()]

rows = {}  # chrom -> [length, bases, min, max]
for p in paths:
    with open(p) as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom = parts[0]
            if chrom == "total":          # drop redundant per-file total row
                continue
            length = int(parts[1]); bases = int(parts[2])
            mn = int(parts[4]); mx = int(parts[5])
            if chrom in rows:             # same chrom seen twice: sum it
                r = rows[chrom]
                r[0] += length; r[1] += bases
                r[2] = min(r[2], mn); r[3] = max(r[3], mx)
            else:
                rows[chrom] = [length, bases, mn, mx]

if not rows:
    raise SystemExit("ERROR: no chromosome rows found in inputs")

def chrom_key(c):
    name = c[3:] if c.lower().startswith("chr") else c
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if name.isdigit():
        return (0, int(name), "")
    if name.upper() in special:
        return (1, special[name.upper()], "")
    return (2, 0, name)

chroms = sorted(rows, key=chrom_key)
tot_len   = sum(rows[c][0] for c in chroms)
tot_bases = sum(rows[c][1] for c in chroms)
tot_mean  = tot_bases / tot_len if tot_len else 0.0
tot_min   = min(rows[c][2] for c in chroms)
tot_max   = max(rows[c][3] for c in chroms)

out_file = sample_id + ".coverage.mosdepth.summary.txt"
with open(out_file, "w") as out:
    out.write("chrom\tlength\tbases\tmean\tmin\tmax\n")
    for c in chroms:
        length, bases, mn, mx = rows[c]
        mean = bases / length if length else 0.0
        out.write(f"{c}\t{length}\t{bases}\t{mean:.2f}\t{mn}\t{mx}\n")
    out.write(f"total\t{tot_len}\t{tot_bases}\t{tot_mean:.2f}\t{tot_min}\t{tot_max}\n")

with open("mean.txt", "w") as m:
    m.write(f"{tot_mean}\n")

print(f"sample={sample_id} chroms={len(chroms)} total_length={tot_len} "
      f"total_bases={tot_bases} mean={tot_mean:.4f}")
PY
  >>>

  output {
    File genome_summary = sample_id + ".coverage.mosdepth.summary.txt"
    Float mean = read_float("mean.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: mem_gb + " GiB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible_tries
    maxRetries: max_retries
  }
}
