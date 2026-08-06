# RNA-seq transcriptome analysis

## Standard order
1. `01_prepare_counts_CPM_PCA.R`: import or load counts and metadata, construct edgeR objects, calculate TMM-normalized CPM values, filter low-expression genes, and perform PCA.
2. `02_edgeR_DEG_GSEA*.R`: run edgeR quasi-likelihood differential-expression analysis and downstream enrichment/GSEA analyses.
3. `03_heatmap*.R`: create heatmaps from expression values and differential-expression gene lists.

## metadata
1. the detail condition of samples can be found in the metadata.rds or metadata.csv