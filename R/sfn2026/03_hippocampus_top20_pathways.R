
# SfN 2026 — Hippocampus top 20 GO pathways
# Uses existing GSEA CSV (no GSEA rerun)

library(dplyr)
library(ggplot2)

# 1. Load existing GSEA results
pathways <- read.csv(
  "results/hippocampus/clusterProfiler_pathways.csv",
  check.names = FALSE
)

# 2. Keep significant pathways
pathways <- pathways %>%
  dplyr::filter(
    !is.na(p.adjust),
    p.adjust < 0.05,
    !is.na(NES)
  )

# 3. Select top 20 unique pathways
top20_terms <- pathways %>%
  dplyr::group_by(Description) %>%
  dplyr::summarise(
    min_padj = min(p.adjust),
    .groups = "drop"
  ) %>%
  dplyr::arrange(min_padj) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::pull(Description)

pathways <- pathways %>%
  dplyr::filter(Description %in% top20_terms)

# 4. Enrichment direction
pathways$Direction <- ifelse(
  pathways$NES > 0,
  "activated",
  "suppressed"
)

# 5. Same comparison order as original plots
comparison_order <- c(
  "Male Dp16 vs WT",
  "Male APP2N vs WT",
  "Male Dp16 vs APP2N",
  "Female Dp16 vs WT",
  "Female APP2N vs WT",
  "Female Dp16 vs APP2N"
)

# Handle original contrast labels if present
comparison_map <- c(
  "Male WTvsDP16" = "Male Dp16 vs WT",
  "Male WTvs2N" = "Male APP2N vs WT",
  "Male DP16v2N" = "Male Dp16 vs APP2N",
  "Female WTvsDP16" = "Female Dp16 vs WT",
  "Female WTvs2N" = "Female APP2N vs WT",
  "Female DP16v2N" = "Female Dp16 vs APP2N"
)

pathways$Comparison <- as.character(pathways$Cluster)

matched <- pathways$Comparison %in% names(comparison_map)

pathways$Comparison[matched] <-
  unname(comparison_map[pathways$Comparison[matched]])

pathways$Comparison <- factor(
  pathways$Comparison,
  levels = unique(c(
    comparison_order,
    pathways$Comparison
  ))
)

# 6. Order and wrap pathway names
pathways$Description <- factor(
  pathways$Description,
  levels = rev(top20_terms),
  labels = vapply(
    rev(top20_terms),
    function(x) paste(strwrap(x, width = 43), collapse = "\n"),
    character(1)
  )
)

# 7. Dot size — use same |NES| fallback
# as the previous CSV figure
pathways$PlotSize <- abs(pathways$NES)

# 8. Plot using original visual formatting
hip_top20 <- ggplot(
  pathways,
  aes(
    x = Comparison,
    y = Description,
    color = p.adjust,
    size = PlotSize
  )
) +
  geom_point() +
  facet_grid(. ~ Direction) +
  scale_color_gradient(
    low = "#e26666",
    high = "#347db5",
    trans = "log10",
    name = "p.adjust"
  ) +
  scale_size_continuous(
    name = "|NES|",
    range = c(3, 8)
  ) +
  labs(
    title = "Hippocampus — Top 20 Significant GO Pathways",
    x = "Comparison",
    y = NULL
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.y = element_text(
      size = 13,
      face = "bold",
      color = "black"
    ),
    axis.text.x = element_text(
      size = 13,
      face = "bold",
      angle = 35,
      hjust = 1,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 18,
      face = "bold",
      margin = margin(t = 15)
    ),
    strip.text.x = element_text(
      size = 22,
      face = "bold",
      color = "black"
    ),
    strip.background = element_rect(
      fill = "grey85",
      color = "black",
      linewidth = 0.8
    ),
    legend.title = element_text(
      size = 16,
      face = "bold"
    ),
    legend.text = element_text(size = 14),
    panel.spacing.x = grid::unit(1.2, "cm"),
    plot.title = element_text(
      size = 22,
      face = "bold"
    ),
    plot.margin = margin(15, 20, 35, 20)
  )

# 9. Preview
print(hip_top20)

# 10. Export
dir.create(
  "figures/sfn2026",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "figures/sfn2026/APPDp16_hippo_top20.pdf",
  plot = hip_top20,
  width = 18,
  height = 11
)

ggsave(
  "figures/sfn2026/APPDp16_hippo_top20.png",
  plot = hip_top20,
  width = 18,
  height = 11,
  dpi = 300,
  bg = "white"
)

cat("Top 20 hippocampus pathway figure saved!\n")
