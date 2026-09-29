library(DESeq2)
library(tidyverse)
library(readxl)
library(tibble)

APPDp16hip <- read.delim("data/hippocampus/APPDp16Hip_counts.txt",
 header = TRUE)

names(APPDp16hip) [1] <- "GeneID"
row.names(APPDp16hip) <- APPDp16hip$GeneID
APPDp16hip <- APPDp16hip [-1]

coldata_raw <- read_xlsx("data/hippocampus/coldata_input.xlsx")

coldata <- as.data.frame(coldata_raw) %>%
  column_to_rownames(var = "Sample") 

coldata$Genotype <- factor(coldata$Genotype)
coldata$Sex <- factor(coldata$Sex)

#find the samples that perfectly match in BOTH your counts and your metadata (needed for DEseq)
common_samples <- intersect(colnames(APPDp16hip), rownames(coldata))

#if this returns 0, sample names are completely different and you will need to manually rename them in your Excel file or R.
length(common_samples) 

# subset both datasets to only include the matching samples,

# ensuring they are in the exact same alphabetical order
APPDp16hip_aligned <- APPDp16hip[, common_samples]
coldata_aligned <- coldata[common_samples, ]

#DESeq2
dds <- DESeqDataSetFromMatrix(countData = APPDp16hip_aligned,
                              colData = coldata_aligned,
                              design = ~ Sex + Genotype + Sex:Genotype)

dds$Genotype <- relevel(dds$Genotype, ref = "WT")

dds <- DESeq(dds)

#for gene lists 
resultsNames(dds)
