library(vsn)
library(pheatmap)

#PCA 
vsd_hyp <- vst(dds_hyp, blind=FALSE)
plotPCA(vsd_hyp, intgroup=c("Sex", "Genotype"))

#heatmap
ntd_hyp <- normTransform(dds_hyp)
meanSdPlot(assay(ntd_hyp))

dds_hyp <- estimateSizeFactors(dds_hyp)
select <- order(rowMeans(counts(dds_hyp, normalized = TRUE)),
                decreasing = TRUE)[1:20]

sample_order_hyp <- order(colData(dds_hyp)$Sex, colData(dds_hyp)$Genotype)
df_hyp <- as.data.frame(colData(dds_hyp)[sample_order_hyp, c("Sex", "Genotype")])

pheatmap(assay(ntd_hyp)[select, sample_order_hyp],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = df_hyp)


#show gene symbols of top 20 

gene_symbols_hyp <- mapIds(edb, keys = rownames(assay(ntd_hyp))[select],
                        column = "SYMBOL", keytype = "GENEID")

pheatmap(assay(ntd_hyp)[select, sample_order_hyp],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         labels_row = gene_symbols_hyp,
         show_colnames = TRUE, #see samples 
         annotation_col = df_hyp)