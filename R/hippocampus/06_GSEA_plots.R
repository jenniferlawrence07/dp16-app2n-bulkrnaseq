library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)

figures_dir <- "figures/hippocampus"
results_dir <- "results/hippocampus"

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ranked gene lists for all 6 comparisons
# Labels reflect the TRUE DESeq2 contrast direction
compare_dp <- list(
  "Male Dp16 vs WT" = geneList_M_WTvsDP16,
  "Male APP2N vs WT" = geneList_M_WTvs2N,
  "Male Dp16 vs APP2N" = geneList_M_DP16vs2N,
  "Female Dp16 vs WT" = geneList_F_WTvsDP16,
  "Female APP2N vs WT" = geneList_F_WTvs2N,
  "Female Dp16 vs APP2N" = geneList_F_DP16vs2N
)

# curated terms for hippocampus
terms_to_keep <- c(
  "cellular response to interferon-beta",
  "myelin sheath",
  "GABA-ergic synapse",
  "glutamate receptor activity",
  "oxidative phosphorylation"
)

# GO GSEA across all 6 comparisons
dpcompare_G_o <- compareCluster(
  geneClusters = compare_dp,
  fun = gseGO,
  ont = "ALL",
  keyType = "ENSEMBL",
  OrgDb = "org.Mm.eg.db",
  verbose = TRUE,
  minGSSize = 3,
  maxGSSize = 800
)

# collapse redundant GO terms
dpcompare_simplified <- clusterProfiler::simplify(
  dpcompare_G_o,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

# keep only curated terms
dpcompare_GO_filtered <- dpcompare_simplified

dpcompare_GO_filtered@compareClusterResult <-
  dpcompare_simplified@compareClusterResult %>%
  dplyr::filter(Description %in% terms_to_keep)

# make plot
APPDp16_hippo_dotplot <-
  dotplot(
    dpcompare_GO_filtered,
    font.size = 16,
    showCategory = 5,
    split = ".sign"
  ) +
  facet_grid(. ~ .sign) +
  labs(
    x = "Comparison",
    y = NULL
  ) +
  theme_bw(base_size = 18) +
  theme(
    # pathway names
    axis.text.y = element_text(
      size = 18,
      face = "bold",
      color = "black"
    ),

    # comparison labels at bottom
    axis.text.x = element_text(
      size = 16,
      face = "bold",
      angle = 35,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),

    # x-axis title
    axis.title.x = element_text(
      size = 18,
      face = "bold",
      margin = margin(t = 15),
      color = "black"
    ),

    # Activated / Suppressed labels
    strip.text.x = element_text(
      size = 22,
      face = "bold",
      color = "black"
    ),

    # legend text
    legend.title = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),

    legend.text = element_text(
      size = 14,
      color = "black"
    ),

    # more breathing room
    panel.spacing.x = unit(1.2, "cm"),

    # margins so labels don't get clipped
    plot.margin = margin(
      t = 15,
      r = 20,
      b = 35,
      l = 20
    )
  )

APPDp16_hippo_dotplot


ggsave(
  file.path(figures_dir, "APPDp16_hippo_dotplot.pdf"),
  plot = APPDp16_hippo_dotplot,
  width = 14,
  height = 8,
  units = "in",
  device = cairo_pdf
)


ggsave(
  file.path(figures_dir, "APPDp16_hippo_dotplot.png"),
  plot = APPDp16_hippo_dotplot,
  width = 14,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)