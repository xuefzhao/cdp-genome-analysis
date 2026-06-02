version 1.0

##
## InferIntervalList.wdl
##
## Inspects a GVCF header to determine which interval list (if any) was
## passed to HaplotypeCaller via -L / --intervals.
##
## Primary method : parse ##GATKCommandLine header lines (HaplotypeCaller,
##                  CombineGVCFs, GenomicsDBImport — checked in that order).
## Fallback       : list chromosomes with data via the tabix index; if all
##                  canonical chroms are present → infer WGS; if only a
##                  subset → infer WES / targeted panel.
##
## Outputs
##   interval_paths  — comma-separated list of -L values found, or "none"
##   interval_type   — "file" | "genomic_regions" | "none (inferred …)"
##   gatk_tool       — GATK tool whose CommandLine was successfully parsed
##   reference       — ##reference value from the VCF header
##   covered_chroms  — chromosomes present in the GVCF tabix index
##   report_txt      — human-readable plain-text report
##

workflow InferIntervalList {

  input {
    File   gvcf
    File   gvcf_index

    String docker    = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    Int    cpu       = 1
    Int    memory_gb = 4
    Int    disk_gb   = 30
  }

  call InferFromHeader {
    input:
      gvcf       = gvcf,
      gvcf_index = gvcf_index,
      docker     = docker,
      cpu        = cpu,
      memory_gb  = memory_gb,
      disk_gb    = disk_gb
  }

  output {
    String interval_paths = InferFromHeader.interval_paths
    String interval_type  = InferFromHeader.interval_type
    String gatk_tool      = InferFromHeader.gatk_tool
    String reference      = InferFromHeader.reference
    String covered_chroms = InferFromHeader.covered_chroms
    File   report_txt     = InferFromHeader.report_txt
  }
}


# ────────────────────────────────────────────────────────────────────────────
task InferFromHeader {
# ────────────────────────────────────────────────────────────────────────────
  input {
    File   gvcf
    File   gvcf_index
    String docker
    Int    cpu
    Int    memory_gb
    Int    disk_gb
  }

  command <<<
    set -euo pipefail
    GVCF="~{gvcf}"

    # ── 1. Extract VCF header ────────────────────────────────────────────────
    bcftools view -h "$GVCF" > header.txt

    # ── 2. List chromosomes that have records (reads index only — fast) ───────
    tabix -l "$GVCF" | tr '\n' ',' | sed 's/,$//' > covered_chroms.txt
    COVERED=$(cat covered_chroms.txt)
    N_CHROMS=$(tabix -l "$GVCF" | wc -l | tr -d ' ')

    # ── 3. Parse header and infer interval list ──────────────────────────────
    # Pass WDL-substituted gvcf path as argv[1]; bash vars as argv[2..3].
    # The heredoc delimiter is quoted ('PYEOF') to suppress bash expansion;
    # WDL's ~{} substitution has already occurred before bash runs.
    python3 - "~{gvcf}" "$COVERED" "$N_CHROMS" << 'PYEOF'
import re, sys

gvcf_path = sys.argv[1]
covered   = sys.argv[2]
n_chroms  = int(sys.argv[3])

with open('header.txt') as fh:
    header = fh.read()

# ── Reference genome ─────────────────────────────────────────────────────────
ref = "not recorded"
for line in header.splitlines():
    if line.startswith('##reference='):
        ref = line.split('=', 1)[1].strip()
        break

# ── Parse GATK CommandLine for -L / --intervals ──────────────────────────────
# GATK 4.x writes:
#   ##GATKCommandLine=<ID=ToolName,CommandLine="ToolName --arg val ...",Version=...>
# We check tools in priority order; stop at the first one found.
TOOLS = [
    'HaplotypeCaller',
    'CombineGVCFs',
    'GenomicsDBImport',
    'SplitIntervals',
    'GenotypeGVCFs',
]

intervals  = []
found_tool = "not found"

for tool in TOOLS:
    pat = r'##GATKCommandLine=<ID=' + re.escape(tool) + r',CommandLine="([^"]+)"'
    m   = re.search(pat, header)
    if not m:
        continue
    found_tool = tool
    cmd = m.group(1)
    # Collect every unique -L / --intervals value (may repeat)
    for iv in re.findall(r'(?:--intervals|-L)\s+(\S+)', cmd):
        if iv not in intervals:
            intervals.append(iv)
    break   # stop at first matched tool

# ── Classify interval type ────────────────────────────────────────────────────
FILE_EXTS = ('.bed', '.interval_list', '.list', '.txt', '.intervals')

def classify(ivs):
    if not ivs:
        return "none"
    if any(
        iv.lower().endswith(FILE_EXTS)
        or iv.startswith('gs://')
        or iv.startswith('s3://')
        or iv.startswith('/')
        for iv in ivs
    ):
        return "file"
    if all(re.match(r'^(chr)?[\dXYMTUW]+(?::\d+-\d+)?$', iv) for iv in ivs):
        return "genomic_regions"
    return "other"

itype = classify(intervals)

# ── When no interval found, infer WGS vs targeted from chrom coverage ─────────
if itype == "none":
    CANONICAL = (
        {f'chr{i}' for i in range(1, 23)} | {'chrX', 'chrY'}
        | {str(i) for i in range(1, 23)} | {'X', 'Y'}
    )
    cov_set = set(covered.split(',')) if covered else set()
    n_canon = len(cov_set & CANONICAL)
    if n_canon >= 22:
        itype = "none (inferred WGS — all canonical chroms covered)"
    elif n_canon >= 3:
        itype = "none (inferred WES/targeted — only %d canonical chroms)" % n_canon
    else:
        itype = "none (inferred targeted panel — very few chroms: %d)" % n_canon

# ── Sample ID from #CHROM line ────────────────────────────────────────────────
sample_id = "unknown"
for line in header.splitlines():
    if line.startswith('#CHROM'):
        parts = line.split('\t')
        if len(parts) > 9:
            sample_id = parts[9].strip()
        break

# ── Write output files ────────────────────────────────────────────────────────
with open('interval_paths.txt', 'w') as f:
    f.write(','.join(intervals) if intervals else 'none')

with open('interval_type.txt', 'w') as f:
    f.write(itype)

with open('gatk_tool.txt', 'w') as f:
    f.write(found_tool)

with open('reference.txt', 'w') as f:
    f.write(ref)

# ── Human-readable report ─────────────────────────────────────────────────────
with open('interval_report.txt', 'w') as f:
    def w(s=""):
        print(s, file=f)
    w("=== Interval List Inference Report ===")
    w("GVCF file       : " + gvcf_path)
    w("Sample ID       : " + sample_id)
    w("Reference       : " + ref)
    w()
    w("GATK tool found : " + found_tool)
    w("Interval type   : " + itype)
    w()
    if intervals:
        w("Interval path(s) / region(s) recorded in header:")
        for iv in intervals:
            w("  " + iv)
    else:
        w("No -L / --intervals argument found in any GATK CommandLine header entry.")
        w("(If the GVCF was generated without -L, this is expected for WGS.)")
    w()
    w("Chromosomes with data in index (%d total):" % n_chroms)
    w("  " + covered)
PYEOF

  >>>

  runtime {
    docker: docker
    cpu:    cpu
    memory: "~{memory_gb} GiB"
    disks:  "local-disk ~{disk_gb} HDD"
  }

  output {
    String interval_paths = read_string("interval_paths.txt")
    String interval_type  = read_string("interval_type.txt")
    String gatk_tool      = read_string("gatk_tool.txt")
    String reference      = read_string("reference.txt")
    String covered_chroms = read_string("covered_chroms.txt")
    File   report_txt     = "interval_report.txt"
  }
}
