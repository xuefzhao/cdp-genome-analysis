# 03_CNV_calling

This directory contains scripts established to detect large coding copy number variation (CNV) from CDP whole-exome sequencing (WES) data using GATK-gCNV.

The full GATK-gCNV pipeline is available on Terra:

https://app.terra.bio/#workspaces/gcnv-dev/gCNV_methods/workflows

This repository includes scripts required to clean and harmonize CDP WES data so the GATK-gCNV pipeline can run successfully.

## Files

- `rename_bam_chroms.py`  
  Python script for renaming chromosome/contig labels in BAM-related inputs to match the naming convention required by downstream processing.

- `rename_bam_chroms.wdl`  
  WDL wrapper for running the chromosome-renaming step in workflow form.

- `rename_bam_chroms_inputs.json`  
  Example input JSON for `rename_bam_chroms.wdl`.

- `.gitkeep`  
  Keeps the directory tracked in Git when needed.
