library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)

figures_dir <- "figures/hypothalamus"
results_dir <- "results/hypothalamus"

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ranked gene lists for all 6 comparisons
# Labels reflect the true DESeq2 contrast direction
compare_dp <- list(
  "Male Dp16 vs WT" = geneList_M_WTvsDP16_hyp,
  "Male APP2N vs WT" = geneList_M_WTvs2N_hyp,
  "Male Dp16 vs APP2N" = geneList_M_DP16vs2N_hyp,
  "Female Dp16 vs WT" = geneList_F_WTvsDP16_hyp,
  "Female APP2N vs WT" = geneList_F_WTvs2N_hyp,
  "Female Dp16 vs APP2N" = geneList_F_DP16vs2N_hyp
)

# curated hypothalamic terms for presentation
terms_to_keep_hypo <- c(
  "cellular response to interferon-beta",
  "myeloid leukocyte activation",
  "mRNA processing",
  "translation at presynapse",
  "regulation of circadian sleep/wake cycle"
)

# GO GSEA across all 6 comparisons
dpcompare_G_o_hyp <- compareCluster(
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
dpcompare_simplified_hyp <- clusterProfiler::simplify(
  dpcompare_G_o_hyp,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

# keep only curated terms
dpcompare_hyp_filtered <- dpcompare_simplified_hyp

dpcompare_hyp_filtered@compareClusterResult <-
  dpcompare_simplified_hyp@compareClusterResult %>%
  dplyr::filter(Description %in% terms_to_keep_hypo)

# rename comparisons for readability
dpcompare_hyp_filtered@compareClusterResult$Cluster <-
  recode(
    dpcompare_hyp_filtered@compareClusterResult$Cluster,
    "Female WTvsDP16" = "F: WT vs Dp16",
    "Female WTvs2N" = "F: WT vs APP2N",
    "Female DP16v2N" = "F: Dp16 vs APP2N",
    "Male WTvsDP16" = "M: WT vs Dp16",
    "Male WTvs2N" = "M: WT vs APP2N",
    "Male DP16v2N" = "M: Dp16 vs APP2N"
  )

# make plot
APPDp16_hyp_dotplot <-
  dotplot(
    dpcompare_hyp_filtered,
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
    text = element_text(
      color = "black"
    ),

    # pathway names
    axis.text.y = element_text(
      size = 18,
      face = "bold",
      color = "black"
    ),

    # comparison labels
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
      color = "black",
      margin = margin(t = 15)
    ),

    # activated / suppressed labels
    strip.text.x = element_text(
      size = 22,
      face = "bold",
      color = "black"
    ),

    # keep gray background behind activated/suppressed
    strip.background = element_rect(
      fill = "grey85",
      color = "black",
      linewidth = 0.8
    ),

    # legend title
    legend.title = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),

    # legend text
    legend.text = element_text(
      size = 14,
      color = "black"
    ),

    # spacing between activated/suppressed panels
    panel.spacing.x = unit(1.2, "cm"),

    # margins so labels don't get clipped
    plot.margin = margin(
      t = 15,
      r = 20,
      b = 35,
      l = 20
    )
  )

# view plot
APPDp16_hyp_dotplot

# high-quality PDF export
ggsave(
  file.path(figures_dir, "APPDp16_hyp_dotplot.pdf"),
  plot = APPDp16_hyp_dotplot,
  width = 14,
  height = 8,
  units = "in",
  device = cairo_pdf
)

# high-quality PNG export
ggsave(
  file.path(figures_dir, "APPDp16_hyp_dotplot.png"),
  plot = APPDp16_hyp_dotplot,
  width = 14,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

# export full simplified pathway table
pathway_df_hyp <- as.data.frame(dpcompare_simplified_hyp)

pathway_df_hyp <- pathway_df_hyp[
  order(
    -pathway_df_hyp$NES,
    pathway_df_hyp$p.adjust
  ),
]

write.csv(
  pathway_df_hyp,
  file = file.path(
    results_dir,
    "clusterProfiler_pathways_simplified_hyp.csv"
  ),
  row.names = FALSE
)


