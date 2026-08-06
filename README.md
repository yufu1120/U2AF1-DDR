# U2AF1-DDR RNA-seq analysis code
This repository contains the cohort-specific code used for Zhang et. al. paper related to U2AF1-HDAC8-DNA Damage Response.

The repository intentionally keeps **separate scripts for each experimental cohort**. It is a record of the analysis code used for the study, not a generalized configuration-driven software package.

## Repository structure

| Directory | Contents |
| --- | --- |
| `RNA-seq-splicing-analysis/` | rMATS SLURM jobs and cohort-specific downstream rMATS analysis scripts |
| `RNA-seq-transcriptome-analysis/` | Count processing, CPM/PCA, edgeR differential-expression/GSEA, and heatmap scripts |
| `environment/` | R, Python, and external-software requirements identified from the scripts |

## Cohorts represented

### Splicing

- MV4-11 and K562, overexpression (HDAC8-OE) and control (EV); ionizing radiation (IR) versus non-irradiated (NIR)
- LSK HDAC8 overexpression (HDAC8-OE) and MIG control; ionizing radiation (IR) versus non-irradiated (NIR)
- LSK Cbfb::MYH11 (CM) and wild-type control; ionizing radiation (IR) versus non-irradiated (NIR)

### Transcriptome

- MV4-11 and K562
- LSK HDAC8-OE
- LSK CM


## General execution order

For splicing analyses:

1. Run the cohort-specific rMATS SLURM jobs in `01-rMATS/`.
2. Run `01_rmats_integrated_pipeline.R` in the corresponding `02-downstream-analysis/` directory.
3. Optionally run the manifest-driven sashimi and relative-frequency-logo utilities.

For transcriptome analyses:

1. Run `01_prepare_counts_CPM_PCA.R`.
2. Run the cohort's `02_edgeR_DEG_GSEA*.R` script.
3. Run `03_heatmap*.R` after the required differential-expression tables exist.

## data source

-Raw and processed RNA-seq data generated in this study are available in the Gene Expression Omnibus (GEO).
-For TARGET AML and Beat AML 1.0, TPM (transcript per million) was derived from the MBatch Omic Browser, the standardized data depository of MD Anderson Cancer Center (https://bioinformatics.mdanderson.org/MQA/), and the clinical information were obtained from the genomic data commons (GDC) portal (https://portal.gdc.cancer.gov/projects/TARGET-AML), and the vizome interactive portal (https://biodev.github.io/BeatAML2/).
-The DDR pathway signature gene list (GO:00006974, 789 genes for human, 772 genes for mouse) used for this analysis were derived from AmiGO 2 (https://amigo.geneontology.org/amigo).

Please contact Dr. Ya-Huei Kuo (ykuo@coh.org) with any questions.
