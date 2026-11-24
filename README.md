# Dp16-App2N Bulk RNA-seq Analysis

**Transcriptomic Analysis of APP Dosage Correction in an Aged Down Syndrome Mouse Model**

This repository contains computational workflows and downstream bulk RNA-seq analyses investigating whether restoring amyloid precursor protein (APP) to wildtype dosage rescues molecular abnormalities in the aged Dp16 mouse model of Down Syndrome.

The project integrates high-performance computing (HPC), differential gene expression analysis, and gene set enrichment analysis to characterize APP-dependent and APP-independent molecular signatures in the hippocampus and hypothalamus.

## 1. Background and Research Question

Down Syndrome (DS), caused by an extra copy of human chromosome 21, is associated with neurological abnormalities, altered sleep architecture, and increased susceptibility to Alzheimer's disease.

APP is triplicated in Down Syndrome and is strongly associated with Alzheimer's disease neuropathology and sleep stability.

Previous studies show that another triplicated gene in DS, RCAN1, partially rescues sleep disturbances in the Dp16 mouse model when it was restored to two copies. 

**Research question:** Does normalizing APP gene dosage rescue transcriptomic abnormalities and alterations in sleep architecture in the hippocampus and hypothalamus of aged Dp16 mice?

We hypothesized that normalizing APP to wildytype levels would restore sleep disturbances and gene programs in aged Dp16 mice. 

## 2. Experimental Design

Bulk RNA sequencing was performed using hippocampal and hypothalamic tissue from aged mice in each of the three genotypes.

| Variable | Description |
| Model | Dp16 mouse model of Down Syndrome (duplication of the region syntenic to Hsa21 on Mmu16)|
| Genotypes | Wild-type (WT), Dp16 (3 copies of APP), Dp16-App2N (2 copies of APP) |
| Age | 10–14 months |
| Sex | Male and female |
| Sample size | n = 4 per sex/genotype |
| Brain regions | Hippocampus and hypothalamus |
| Sequencing approach | Bulk RNA-seq |
| Differential expression | DESeq2 |
| Pathway analysis | Gene Set Enrichment Analysis (GSEA) |

Analyses were performed independently by brain region, with comparisons across genotype and sex.

### Primary Comparisons

1. **WT vs. Dp16:** Identify transcriptomic abnormalities associated with trisomy.
2. **WT vs. Dp16-App2N:** Determine which abnormalities remain following APP dosage normalization.
3. **Dp16 vs. Dp16-App2N:** Identify molecular differences associated with restoring APP to disomic dosage.

## 3. Computational Workflow

The analysis consists of upstream sequencing processing on CU Boulder's Alpine HPC cluster and downstream transcriptomic analysis in R.

### Upstream RNA-seq Processing

Raw sequencing reads were processed using the following tools:

1. **fastp:** Paired-end read quality control and trimming.
2. **STAR:** Alignment to the mouse reference genome. (mm39)
3. **featureCounts:** Gene-level read quantification.

Processing scripts are organized by brain region and use SLURM for HPC job submission.

### Downstream Transcriptomic Analysis

Downstream analyses were performed in R and include:

- Count matrix preparation and QC
- Principal component analysis (PCA)
- Sample clustering and heatmap visualization
- Differential gene expression using DESeq2
- Ranked gene list generation (DESeq2 Wald stat)
- Gene set enrichment analysis using clusterProfiler
- Volcano plots and enrichment visualizations

**Statistical thresholds:**
- Differentially expressed genes: adjusted p-value < 0.1
- Gene set enrichment: FDR < 0.05

## 4. Key Findings

### Differential Gene Expression

Comparisons between WT and Dp16 identified approximately **97–120 differentially expressed genes (DEGs)** across the analyzed sex and brain-region groups.

WT versus Dp16-App2N comparisons identified **164–340 DEGs**, whereas direct comparisons between Dp16 and Dp16-App2N identified only **0–1 DEG per group**.

These findings indicate limited detectable normalization of the transcriptome following APP dosage correction.

### Persistent Molecular Abnormalities

Gene set enrichment analysis identified biological pathways that remained altered in both Dp16 and Dp16-App2N mice relative to WT.

Persistent molecular signatures included:

- Immune and interferon-associated processes
- Calcium-associated signaling
- Suppression of hippocampal RNA-processing programs
- Reduced hypothalamic chromatin-remodeling pathways
- Reduced histone acetyltransferase-complex pathways

Response to interferon-beta remained strongly enriched across brain regions and sexes, with normalized enrichment scores of 2.05–2.40 and adjusted p-values below 0.0001 in the strongest comparisons.

### Interpretation

Restoring APP to disomic dosage did not broadly normalize gene expression or pathway-level abnormalities in aged Dp16 mice.

Together with persistent sleep abnormalities, these exploratory findings suggest that APP-independent molecular mechanisms may contribute to neurological dysfunction associated with Down Syndrome.

## 5. Repository Structure

```text
dp16-app2n-bulkrnaseq/
│
├── HPC/
│   ├── hippocampus/
│   │   ├── fastp_hippo.sh
│   │   ├── star_hippo.sh
│   │   ├── featureCounts.R
│   │   └── submit_feature_counts.sh
│   │
│   └── hypothalamus/
│       ├── fastp_hyp.sh
│       ├── star_hyp.sh
│       ├── featureCounts.R
│       └── submit_feature_counts.sh
│
├── R/
│   └── Region-specific downstream analysis scripts
│
├── data/
│   └── Processed count matrices and metadata
│
├── results/
│   └── Differential expression and enrichment results
│
├── figures/
│   └── Generated analysis figures
│
├── docs/
│   └── Supporting documentation
│
└── README.md
```

## 6. Software and Computational Tools

| Category | Tools |
| HPC environment | CU Boulder Alpine, Linux, SLURM |
| Quality control and trimming | fastp |
| Read alignment | STAR |
| Gene quantification | featureCounts |
| Statistical analysis | R, DESeq2 |
| Functional enrichment | clusterProfiler, GSEA |
| Visualization | R-based PCA, heatmaps, volcano plots, enrichment figures |
| Version control | Git, GitHub |

## 7. Reproducibility

The repository is organized to separate sequencing processing, downstream analyses, and generated results.

### HPC Processing

The `HPC/` directory contains brain-region-specific processing scripts originally designed for execution on CU Boulder's Alpine cluster.

Processing follows the general sequence:

**FASTQ → fastp → STAR → featureCounts → Gene Count Matrix**

The scripts include HPC-specific reference paths, input/output directories, and SLURM settings that must be configured for the target computing environment.

The HPC scripts are provided as research processing code and have not yet been validated as a fully automated, end-to-end workflow.

### Downstream Analysis

The `R/` directory contains scripts for differential expression and pathway analysis.

Processed count matrices and associated sample metadata are organized in `data/`, with analytical outputs in `results/` and `figures/`.

Relative file paths are used where possible to support execution from the repository root.

Dependencies, sample metadata, and analysis-specific inputs should be checked before reproducing individual results.

### Ongoing Development

Future improvements include:

- Standardizing software environments and dependency management
- Implementing a modular Snakemake workflow for HPC processing
- Automating sample handling and workflow dependencies
- Improving reproducibility through testing and documentation

## 8. Research Context

This repository supports exploratory research examining the molecular consequences of APP gene dosage correction in an aged mouse model of Down Syndrome.

**Associated scientific abstract:**

Lawrence, J., Milstead, R., Cain, P., Borski, C., Zou, S., Kastengren, K., Kushner, K., & Hoeffer, C. (2026).Restoring APP to disomic dosage does not rescue sleep and reveals persistent APP-independent molecular abnormalities in aged Dp16 mice. Society for Neuroscience (SfN) Annual Meeting, 2026.

## 9. Project Scope

This repository documents computational methods and exploratory transcriptomic findings associated with an ongoing neuroscience research project.

The analyses are intended to support reproducible research and investigation of gene-regulatory mechanisms in Down Syndrome. The findings should not be interpreted as establishing causal effects of individual enriched pathways without experimental validation.
