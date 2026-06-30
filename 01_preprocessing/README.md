# 01_preprocessing

Preprocessing and quality-control scripts for the **cdp-genome-analysis** project,
focused on whole-exome sequencing (WES) samples from individuals affected by
**constitutional delay of puberty (CDP)**, collected from a variety of sources.

## Purpose

The scripts in this folder cover two preprocessing steps:

1. **Library / coverage quality evaluation** — assess the sequencing library
   quality and coverage of each WES sample. Coverage is computed with mosdepth,
   the per-chromosome outputs are integrated into a genome-wide summary, and
   on-target vs. off-target depth is summarized. Because samples are drawn from
   multiple resources, this provides a consistent QC view across cohorts.
2. **Genetic ancestry inference** — infer the genetic ancestry of each sample.

## Workflows

| Step | Workflow | Description |
|------|----------|-------------|
| Coverage QC | [`MosDepth.wdl`](MosDepth.wdl) | Runs [mosdepth](https://github.com/brentp/mosdepth) to evaluate per-sample coverage of the WES data (run per chromosome). |
| Coverage QC | [`IntegrateMosdepthSummary.wdl`](IntegrateMosdepthSummary.wdl) | Integrates a sample's per-chromosome mosdepth `*.summary.txt` files into one genome-wide summary: sums `length` and `bases` across chromosomes and recomputes the genome **mean** coverage (`sum(bases)/sum(length)`). Outputs the combined summary file and the recalculated mean. |
| Coverage QC | [`TargetDepthSummary.wdl`](TargetDepthSummary.wdl) | Summarizes per-base depth **within** vs. **outside** the WES capture targets (mean / median / sd), using mosdepth per-base coverage and the target intervals BED. |
| Ancestry inference | [`ancestry-inference-hail-v01.wdl`](ancestry-inference-hail-v01.wdl) | Infers genetic ancestry for each sample (Hail-based pipeline). |

## Notes

- `mosdepth` is used to summarize coverage/depth as a measure of WES library quality.
- `IntegrateMosdepthSummary.wdl` recomputes the genome mean as a base-weighted
  value (`sum(bases)/sum(length)`), not a naive average of per-chromosome means,
  and drops the redundant `total` row present in each per-chromosome file.
- Ancestry is inferred per sample via the Hail-based workflow above.
