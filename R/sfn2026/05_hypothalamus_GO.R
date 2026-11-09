
# SfN 2026 — Hypothalamus GO pathways
# Uses saved GSEA results; does NOT rerun GSEA

library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(ggplot2)

figures_dir <- "figures/sfn2026"
results_dir <- "results/sfn2026"

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load existing simplified GSEA results
gsea_file <- file.path(
  results_dir,
  "GSEA_hypothalamus_simplified.rds"
)

if (file.exists(gsea_file)) {
  dpcompare_simplified_hyp <- readRDS(gsea_file)
} else if (!exists("dpcompare_simplified_hyp")) {
  stop("Saved hypothalamus GSEA object not found.")
}

# 2. Eight selected pathways
terms_to_keep_hypo <- c(
  "cellular response to interferon-beta",
  "myeloid leukocyte activation",
  "mRNA processing",
  "translation at presynapse",
  "regulation of circadian sleep/wake cycle",
  "chromatin remodeling",
  "histone acetyltransferase complex",
  "regulation of circadian rhythm"
)

# 3. Filter saved GSEA results
dpcompare_hyp_filtered <- dpcompare_simplified_hyp

dpcompare_hyp_filtered@compareClusterResult <-
  dpcompare_simplified_hyp@compareClusterResult %>%
  dplyr::filter(
    Description %in% terms_to_keep_hypo,
    !is.na(p.adjust),
    p.adjust < 0.05
  )

found_terms <- unique(
  dpcompare_hyp_filtered@compareClusterResult$Description
)

cat("Hypothalamus pathways included:\n")
print(found_terms)

missing_terms <- setdiff(terms_to_keep_hypo, found_terms)

if (length(missing_terms) > 0) {
  warning(
    "Missing pathways: ",
    paste(missing_terms, collapse = ", ")
  )
}

# 4. Short comparison labels
short_labels <- c(
  "Male Dp16 vs WT" = "M Dp16/WT",
  "Male APP2N vs WT" = "M APP2N/WT",
  "Male Dp16 vs APP2N" = "M Dp16/2N",
  "Female Dp16 vs WT" = "F Dp16/WT",
  "Female APP2N vs WT" = "F APP2N/WT",
  "Female Dp16 vs APP2N" = "F Dp16/2N"
)

format_labels <- function(x) {
  x <- as.character(x)
  labels <- unname(short_labels[x])
  labels[is.na(labels)] <- x[is.na(labels)]
  labels
}

# 5. Original clusterProfiler dotplot styling
APPDp16_hyp_dotplot <- dotplot(
  dpcompare_hyp_filtered,
  font.size = 14,
  showCategory = 100,
  split = ".sign"
) +
  facet_grid(. ~ .sign) +
  scale_x_discrete(labels = format_labels) +
  labs(
    x = "Comparison",
    y = NULL
  ) +
  theme_bw(base_size = 16) +
  theme(
    text = element_text(color = "black"),

    axis.text.y = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),

    axis.text.x = element_text(
      size = 10,
      face = "bold",
      angle = 60,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),

    axis.title.x = element_text(
      size = 15,
      face = "bold",
      color = "black",
      margin = margin(t = 10)
    ),

    strip.text.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),

    strip.background = element_rect(
      fill = "grey85",
      color = "black",
      linewidth = 0.8
    ),

    legend.title = element_text(
      size = 11,
      face = "bold",
      color = "black"
    ),

    legend.text = element_text(
      size = 10,
      color = "black"
    ),

    panel.spacing.x = grid::unit(0.5, "cm"),

    plot.margin = margin(
      t = 5,
      r = 5,
      b = 10,
      l = 5
    )
  )

# 6. Preview
print(APPDp16_hyp_dotplot)

# Verify plotted pathways
cat("\nPathways displayed in hypothalamus plot:\n")
print(unique(as.character(
  APPDp16_hyp_dotplot$data$Description
)))

# 7. Save PDF
ggsave(
  file.path(figures_dir, "APPDp16_hyp_GO_8terms.pdf"),
  plot = APPDp16_hyp_dotplot,
  width = 10,
  height = 9,
  units = "in",
  device = "pdf"
)

# 8. Save PNG
ggsave(
  file.path(figures_dir, "APPDp16_hyp_GO_8terms.png"),
  plot = APPDp16_hyp_dotplot,
  width = 10,
  height = 9,
  units = "in",
  dpi = 600,
  bg = "white"
)

cat("\nHypothalamus GO figure saved successfully!\n")
