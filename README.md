# Dp16-App2N Bulk RNA-seq

Bulk RNA-seq analysis for the Dp16 / APP dosage correction project.

## Brain regions
- Hippocampus
- Hypothalamus

## Main analyses
- DESeq2
- PCA / heatmaps
- Ranked gene lists
- Volcano plots
- GSEA / clusterProfiler

## Project structure
- `R/` analysis scripts
- `HPC/` reserved for Alpine / sequencing-processing scripts; this part of the workflow is not currently developed/actively used
- `data/` count matrices and metadata
- `results/` exported analysis results
- `figures/` generated figures
- `docs/` supporting project documents

Relative file paths are used so the same code can run from different computers when the repository root is the working directory.
