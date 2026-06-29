# 01_preprocessing

Preprocessing and quality-control scripts for the **cdp-genome-analysis** project,
focused on whole-exome sequencing (WES) samples from individuals affected by
**constitutional delay of puberty (CDP)**, collected from a variety of sources.

## Purpose

The scripts in this folder cover two preprocessing steps:

1. **Library / coverage quality evaluation** — assess the sequencing library
   quality and coverage of each WES sample. Samples are drawn from multiple
   resources, so this step provides a consistent QC view across cohorts.
2. **Genetic ancestry inference** — infer the genetic ancestry of each sample.

## Workflows

| Step | Workflow | Description |
|------|----------|-------------|
| Coverage QC | [`MosDepth.wdl`](MosDepth.wdl) | Runs [mosdepth](https://github.com/brentp/mosdepth) to evaluate per-sample coverage of the WES data. |
| Ancestry inference | [`ancestry-inference-hail-v01.wdl`](ancestry-inference-hail-v01.wdl) | Infers genetic ancestry for each sample (Hail-based pipeline). |

## Notes

- `mosdepth` is used to summarize coverage/depth as a measure of WES library quality.
- Ancestry is inferred per sample via the Hail-based workflow above.
