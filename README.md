# CDP Genome Analysis Pipeline

Scripts and workflows for processing whole-genome sequencing (WGS) data from
**Constitutional Delayed Puberty (CDP)** samples — from raw sequencing reads
through joint-called variant files and downstream statistical analyses.

---

## Overview

This repository contains the complete bioinformatics pipeline used to:

1. Process raw WGS data (FASTQ → aligned BAM)
2. Call structural variants (SVs), short tandem repeats (STRs), and small
   variants (SNVs/indels) per sample
3. Jointly genotype variants across all CDP samples
4. Annotate variants with functional consequences and population allele frequencies
5. Perform statistical analyses to identify variants and genes associated
   with delayed puberty phenotypes

---

## Repository Structure

```
cdp-genome-analysis/
├── 01_preprocessing/          # Raw data QC, alignment, and BAM processing
├── 02_variant_calling/        # Per-sample variant calling (SVs, STRs, SNVs)
├── 03_joint_calling/          # Joint genotyping across all CDP samples
├── 04_annotation/             # Variant annotation (VEP, gnomAD AF, etc.)
├── 05_statistical_analysis/   # Burden tests, association analyses, plots
└── utils/                     # Shared helper scripts and utilities
```

---

## Pipeline Stages

### 01 — Preprocessing
- FASTQ quality control (FastQC / MultiQC)
- Read alignment to GRCh38 (BWA-MEM2 / DRAGEN)
- BAM sorting, duplicate marking, and base quality score recalibration (BQSR)
- Coverage and alignment QC metrics

### 02 — Variant Calling (per sample)
- **SVs**: GATK-SV, Manta, PBSV (long reads), PBMM2
- **STRs / tandem repeats**: ExpansionHunter, TRGT, HIPHASE
- **SNVs/indels**: DeepVariant, GATK HaplotypeCaller
- Long-read assembly and PAV for phased variant detection

### 03 — Joint Calling
- SV joint genotyping across all CDP samples (GATK-SV module 04–07)
- SNV/indel joint genotyping via GATK GenotypeGVCFs
- STR joint calling and phasing
- Hard filtering and VQSR for small variant cohort callset

### 04 — Annotation
- Functional annotation with VEP (Ensembl)
- Population allele frequencies from gnomAD v4 (short reads) and gnomAD LR
- CADD, SpliceAI, pLI/LOEUF scores
- Overlap with known puberty/HPG-axis loci and OMIM disease genes

### 05 — Statistical Analysis
- Rare variant burden tests (case vs. control)
- Gene-based collapsing analyses (pLoF, missense, regulatory)
- STR/TR expansion burden in CDP vs. general population
- Gene family / segmental duplication burden (paralog-aware)
- Visualization and result reporting

---

## Data

| Category        | Description                              |
|-----------------|------------------------------------------|
| Cohort          | Constitutional Delayed Puberty (CDP)     |
| Sequencing      | Whole-genome sequencing (short + long read) |
| Reference       | GRCh38                                   |
| Sample size     | TBD                                      |

---

## Dependencies

- **Workflow engine**: WDL / Cromwell (Terra cloud platform)
- **Alignment**: BWA-MEM2, PBMM2
- **SV calling**: GATK-SV, Manta, PBSV
- **STR calling**: ExpansionHunter, TRGT
- **SNV calling**: DeepVariant, GATK HaplotypeCaller
- **Annotation**: VEP, bcftools, pysam, pandas
- **Statistics**: R (burden tests), Python (scipy, statsmodels)

---

## Contact

Xuefei Zhao — xuefzhao@umich.edu  
Talkowski Lab, Massachusetts General Hospital / Harvard Medical School
