library(DESeq2)
library(tidyverse)
library(readxl)
library(tibble)

APPDp16hyp<- read.delim("data/hypothalamus/APPDp16Hyp_counts.txt",
 header = TRUE)

names(APPDp16hyp) [1] <- "GeneID"
row.names(APPDp16hyp) <- APPDp16hyp$GeneID
APPDp16hyp <- APPDp16hyp [-1]


coldata_raw_hyp <- read_xlsx("data/hypothalamus/hyp_coldata_input.xlsx")

coldata_hyp <- as.data.frame(coldata_raw_hyp) %>%
  column_to_rownames(var = "Sample") 

coldata_hyp$Genotype <- factor(coldata_hyp$Genotype)
coldata_hyp$Sex <- factor(coldata_hyp$Sex)

#find the samples that perfectly match in BOTH your counts and your metadata (needed for DEseq)
common_samples_hyp <- intersect(colnames(APPDp16hyp), rownames(coldata_hyp))

#if this returns 0, sample names are completely different and you will need to manually rename them in your Excel file or R.
length(common_samples_hyp) 

# subset both datasets to only include the matching samples,

# ensuring they are in the exact same alphabetical order
APPDp16hyp_aligned <- APPDp16hyp[, common_samples_hyp]
coldata_aligned_hyp <- coldata_hyp[common_samples_hyp, ]

#DESeq2
dds_hyp <- DESeqDataSetFromMatrix(countData = APPDp16hyp_aligned,
                              colData = coldata_aligned_hyp,
                              design = ~ Sex + Genotype + Sex:Genotype)

dds_hyp$Genotype <- relevel(dds_hyp$Genotype, ref = "WT")

dds_hyp <- DESeq(dds_hyp)

#for gene lists 
resultsNames(dds_hyp)