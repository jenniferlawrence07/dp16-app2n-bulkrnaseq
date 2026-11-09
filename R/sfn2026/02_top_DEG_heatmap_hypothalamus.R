
# SfN 2026 — Top DEG heatmap
# Hypothalamus
# Matches the hippocampus heatmap design

library(DESeq2)
library(dplyr)
library(pheatmap)
library(AnnotationDbi)
library(EnsDb.Mmusculus.v79)
library(ensembldb)

# 1. Load saved hypothalamus DESeq2 dataset
dds_hyp <- readRDS(
  "results/sfn2026/dds_hypothalamus.rds"
)

# Retrieve existing sex-specific DESeq2 contrasts
# This does NOT refit the DESeq2 model.
source("R/hypothalamus/02_gene_lists.R")

deg_results <- list(
  Male_WT_vs_Dp16 = res_M_WTvsDP16_hyp,
  Male_WT_vs_APP2N = res_M_WTvs2N_hyp,
  Female_WT_vs_Dp16 = res_F_WTvsDP16_hyp,
  Female_WT_vs_APP2N = res_F_WTvs2N_hyp
)

# 2. Select the four WT comparisons
contrasts <- c(
  "Male_WT_vs_Dp16",
  "Male_WT_vs_APP2N",
  "Female_WT_vs_Dp16",
  "Female_WT_vs_APP2N"
)

# 3. Genes flagged in the hippocampus heatmap
# Apply the same visualization exclusions here.
# This does not modify the underlying DESeq2 results.
excluded_gene_ids <- c(
  "ENSMUSG00000076583",  # Igkv8-24
  "ENSMUSG00000094420"   # Igkv10-96
)

# Select top 10 significant genes per comparison
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
  unlist(
    lapply(deg_results[contrasts], get_top_genes),
    use.names = FALSE
  )
)

cat(
  "Number of unique hypothalamic genes:",
  length(selected_genes),
  "\n"
)

if (length(selected_genes) < 2) {
  stop("Fewer than two genes selected for heatmap.")
}

# 4. Variance-stabilizing transformation
vsd <- vst(dds_hyp, blind = FALSE)

expression_matrix <- assay(vsd)[
  selected_genes, ,
  drop = FALSE
]

# 5. Organize mice by sex and genotype
metadata <- as.data.frame(colData(dds_hyp))

metadata$Sex <- factor(
  metadata$Sex,
  levels = c("Female", "Male")
)

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

# Remove any rows with undefined Z-scores
valid_rows <- apply(scaled_matrix, 1, function(x) {
  all(is.finite(x))
})

scaled_matrix <- scaled_matrix[
  valid_rows, ,
  drop = FALSE
]

# 7. Annotate with mouse gene symbols
edb <- EnsDb.Mmusculus.v79

gene_ids <- rownames(scaled_matrix)

gene_symbols <- AnnotationDbi::mapIds(
  x = edb,
  keys = gene_ids,
  column = "SYMBOL",
  keytype = "GENEID",
  multiVals = "first"
)

gene_labels <- as.character(gene_symbols)

# Keep Ensembl IDs for genes without symbols
missing_labels <- is.na(gene_labels) | gene_labels == ""

gene_labels[missing_labels] <- gene_ids[missing_labels]

# Prevent duplicate labels
rownames(scaled_matrix) <- make.unique(gene_labels)

# 8. Identify mouse chromosome 16 genes
chr16_genes <- genes(
  edb,
  filter = AnnotationFilter::SeqNameFilter("16")
)$gene_id

# Annotate rows using original Ensembl IDs
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

# Verify gene annotations
cat("Mapped hypothalamic gene labels:\n")
print(rownames(scaled_matrix))

# 9. Generate and save heatmap
output <- "figures/sfn2026/top_DEGs_hypothalamus"

dir.create(
  "figures/sfn2026",
  recursive = TRUE,
  showWarnings = FALSE
)

# Same colors as hippocampus heatmap
heatmap_colors <- colorRampPalette(
  c("#2600ff", "white", "#ff2301")
)(100)

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
  main = "Hypothalamus — Top Differentially Expressed Genes",
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

cat("Hypothalamus heatmap saved successfully!\n")
