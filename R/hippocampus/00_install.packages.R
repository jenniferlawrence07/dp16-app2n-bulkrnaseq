if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("tidyverse", "pheatmap", "ggrepel", "ggridges", "openxlsx", "ggplot2")
bioc_pkgs <- c(
  "DESeq2", "vsn", "AnnotationDbi", "EnsDb.Mmusculus.v79",
  "clusterProfiler", "pathview", "enrichplot", "org.Mm.eg.db", "DOSE"
)
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
}
message("All packages present.")
