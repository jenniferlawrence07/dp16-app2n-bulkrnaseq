
# SfN 2026 — Top DEG heatmap
# Hippocampus

library(DESeq2)
library(dplyr)
library(pheatmap)
library(AnnotationDbi)
library(EnsDb.Mmusculus.v79)

# 1. Load saved analysis objects
dds <- readRDS("results/sfn2026/dds_hippocampus.rds")

deg_results <- readRDS(
  "results/sfn2026/DEG_results_hippocampus.rds"
)

# 2. Select the four WT comparisons
contrasts <- c(
  "Male_WT_vs_Dp16",
  "Male_WT_vs_APP2N",
  "Female_WT_vs_Dp16",
  "Female_WT_vs_APP2N"
)

# Previously flagged sparse-expression gene
excluded_gene_ids <- c(
  "ENSMUSG00000076583",  # Igkv8-24
  "ENSMUSG00000094420"   # Igkv10-96
)

get_top_genes <- function(res, n = 10) {
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)

  df %>%
    dplyr::filter(
      !is.na(padj),
      padj < 0.1,
      !gene_id %in% excluded_gene_ids
    ) %>%
    dplyr::arrange(
      padj,
      dplyr::desc(abs(log2FoldChange))
    ) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::pull(gene_id)
}


selected_genes <- unique(
  unlist(lapply(
    deg_results[contrasts],
    get_top_genes
  ), use.names = FALSE)
)

cat("Number of unique genes:", length(selected_genes), "\n")

# 4. Variance-stabilizing transformation
vsd <- vst(dds, blind = FALSE)

expression_matrix <- assay(vsd)[selected_genes, , drop = FALSE]

# 5. Organize mice by sex and genotype
metadata <- as.data.frame(colData(dds))

metadata$Sex <- factor(metadata$Sex)
metadata$Genotype <- factor(
  metadata$Genotype,
  levels = c("WT", "Dp16", "APP2N")
)

sample_order <- rownames(metadata)[
  order(metadata$Sex, metadata$Genotype)
]

expression_matrix <- expression_matrix[, sample_order]

annotation_col <- metadata[
  sample_order,
  c("Sex", "Genotype"),
  drop = FALSE
]


# 6. Scale each gene to a row Z-score
scaled_matrix <- t(scale(t(expression_matrix)))

# 7. Annotate with mouse gene symbols
# Use the same annotation database as the original volcano analysis

edb <- EnsDb.Mmusculus.v79

gene_ids <- rownames(scaled_matrix)

gene_symbols <- AnnotationDbi::mapIds(
  x = edb,
  keys = gene_ids,
  column = "SYMBOL",
  keytype = "GENEID",
  multiVals = "first"
)

# Convert to an ordinary character vector
gene_labels <- as.character(gene_symbols)

# Keep Ensembl IDs for genes without an annotation
missing_labels <- is.na(gene_labels) | gene_labels == ""

gene_labels[missing_labels] <- gene_ids[missing_labels]

# Prevent duplicate row names
rownames(scaled_matrix) <- make.unique(gene_labels)

# Identify mouse chromosome 16 genes
chr16_genes <- genes(
  edb,
  filter = AnnotationFilter::SeqNameFilter("16")
)$gene_id

# Annotate heatmap rows using original Ensembl IDs
row_annotation <- data.frame(
  Chromosome = ifelse(
    gene_ids %in% chr16_genes,
    "Mmu16",
    "Other"
  )
)

rownames(row_annotation) <- rownames(scaled_matrix)

# Annotation colors
annotation_colors <- list(
  Chromosome = c(
    Mmu16 = "#e78c04",
    Other = "#B7C0CC"
  )
)

# Verify annotation
cat("Mapped gene labels:\n")
print(rownames(scaled_matrix))


# 8. Generate and save the annotated heatmap

output <- "figures/sfn2026/top_DEGs_hippocampus"

dir.create(
  "figures/sfn2026",
  recursive = TRUE,
  showWarnings = FALSE
)

# Heatmap color palette
heatmap_colors <- colorRampPalette(
  c("#2600ff", "white", "#ff2301")
)(100)

# Generate heatmap once
heatmap_plot <- pheatmap(
  scaled_matrix,
  annotation_col = annotation_col,
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  annotation_names_row = FALSE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  fontsize_row = 10,
  main = "Hippocampus — Top Differentially Expressed Genes",
  color = heatmap_colors,
  breaks = seq(-3, 3, length.out = 101),
  border_color = NA,
  silent = TRUE
)

# Preview in R
grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)

# Save PDF
pdf(paste0(output, ".pdf"), width = 12, height = 9)
grid::grid.draw(heatmap_plot$gtable)
dev.off()

# Save high-resolution PNG
png(
  paste0(output, ".png"),
  width = 3600,
  height = 2700,
  res = 300
)
grid::grid.draw(heatmap_plot$gtable)
dev.off()

cat("Annotated hippocampus heatmap saved successfully!\n")
