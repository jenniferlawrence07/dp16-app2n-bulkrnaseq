library(ggplot2)
library(ggrepel)
library(dplyr)
library(ensembldb)

# pull chr16 gene IDs 
library(ensembldb)
chr16_genes <- genes(edb, filter = AnnotationFilter::SeqNameFilter("16"))$gene_id

PADJ_CUTOFF <- 0.1
BASEMEAN_CUTOFF <- 11

volcano_plot <- function(plot_df, plot_title,
                          padj_cutoff = PADJ_CUTOFF,
                          basemean_cutoff = BASEMEAN_CUTOFF) {

  plot_df <- plot_df %>% dplyr::filter(baseMean >= basemean_cutoff) # nolint: object_usage_linter.

  # define significance
  plot_df$Significance <- "Not Significant"
  plot_df$Significance[plot_df$padj < padj_cutoff &
                         plot_df$log2FoldChange > 0 &
                         !is.na(plot_df$padj)] <- "Up-regulated"
  plot_df$Significance[plot_df$padj < padj_cutoff &
                         plot_df$log2FoldChange < 0 &
                         !is.na(plot_df$padj)] <- "Down-regulated"

  # chr16 flag: label up/down genes that are also on chr16
  plot_df$Color_Group <- plot_df$Significance
  plot_df$Color_Group[plot_df$Significance == "Up-regulated" &
                        plot_df$ensembl_id %in% chr16_genes] <- "Up-regulated (Chr16)"
  plot_df$Color_Group[plot_df$Significance == "Down-regulated" &
                        plot_df$ensembl_id %in% chr16_genes] <- "Down-regulated (Chr16)"

  # top 10 per group for labels
  top_up <- plot_df %>% dplyr::filter(Color_Group == "Up-regulated") %>% arrange(padj) %>% slice_head(n = 10)
  top_up_chr16 <- plot_df %>% dplyr::filter(Color_Group == "Up-regulated (Chr16)") %>% arrange(padj) %>% slice_head(n = 10)
  top_down <- plot_df %>% dplyr::filter(Color_Group == "Down-regulated") %>% arrange(padj) %>% slice_head(n = 10)
  top_down_chr16 <- plot_df %>% dplyr::filter(Color_Group == "Down-regulated (Chr16)") %>% arrange(padj) %>% slice_head(n = 10)

  label_df <- bind_rows(top_up, top_up_chr16, top_down, top_down_chr16) %>%
    mutate(label_color = case_when(
      Color_Group == "Up-regulated (Chr16)" ~ "red4",
      Color_Group == "Up-regulated" ~ "firebrick1",
      Color_Group == "Down-regulated (Chr16)" ~ "royalblue4",
      Color_Group == "Down-regulated" ~ "royalblue1",
      TRUE ~ "grey40"
    ))

  p <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Color_Group)) +
    geom_point(alpha = 0.8, size = 1.5) +
    scale_color_manual(
      name = paste0("Significance (padj < ", padj_cutoff, ")"),
      values = c(
        "Not Significant" = "grey80",
        "Up-regulated" = "firebrick1",
        "Down-regulated" = "royalblue1",
        "Up-regulated (Chr16)" = "red4",
        "Down-regulated (Chr16)" = "royalblue4"
      ),
      guide = guide_legend(override.aes = list(size = 3, alpha = 1))
    ) +
    geom_text_repel(
      data = label_df,
      aes(label = symbol, color = I(label_color)),
      size = 3,
      max.overlaps = 30,
      box.padding = 0.5,
      segment.color = "grey40",
      segment.size = 0.3,
      fontface = "plain",
      show.legend = FALSE
    ) +
    theme_minimal() +
    labs(title = plot_title, x = "Log2 Fold Change", y = "-Log10(p-value)") +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "plain"),
      legend.position = "right",
      legend.title    = element_text(size = 8, face = "plain"),
      legend.text     = element_text(size = 7),
      legend.key.size = unit(0.4, "cm")
    )
  figures_dir <- "figures/hypothalamus"
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  filename <- paste0(gsub("[^A-Za-z0-9]+", "_", plot_title), ".png")

  ggsave(filename = file.path(figures_dir, filename), plot = p,
         dpi = 600, width = 8, height = 6, units = "in")

  return(p)
}

# PIPELINE
## MALES
volcano_plot(volcano_M_WTvsDP16_hyp, "Male WT vs Dp16")
volcano_plot(volcano_M_WTvs2N_hyp,   "Male WT vs APP2N")
volcano_plot(volcano_M_DP16vs2N_hyp, "Male Dp16 vs APP2N")

## FEMALES
volcano_plot(volcano_F_WTvsDP16_hyp, "Female WT vs Dp16")
volcano_plot(volcano_F_WTvs2N_hyp,   "Female WT vs APP2N")
volcano_plot(volcano_F_DP16vs2N_hyp, "Female Dp16 vs APP2N")