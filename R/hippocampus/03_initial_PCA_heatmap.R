
library(vsn)
library(pheatmap)

#PCA 
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("Sex", "Genotype"))


#heatmap
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

dds <- estimateSizeFactors(dds)
select <- order(rowMeans(counts(dds, normalized = TRUE)),
                decreasing = TRUE)[1:20]

sample_order <- order(colData(dds)$Sex, colData(dds)$Genotype)
df <- as.data.frame(colData(dds)[sample_order, c("Sex", "Genotype")])

pheatmap(assay(ntd)[select, sample_order],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = df)


#show gene symbols of top 20 

gene_symbols <- mapIds(edb, keys = rownames(assay(ntd))[select],
                        column = "SYMBOL", keytype = "GENEID")

pheatmap(assay(ntd)[select, sample_order],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         labels_row = gene_symbols,
         show_colnames = TRUE, #see samples 
         annotation_col = df)