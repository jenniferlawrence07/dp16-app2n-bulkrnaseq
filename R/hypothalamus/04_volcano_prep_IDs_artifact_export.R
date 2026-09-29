library(dplyr)
library(AnnotationDbi)
library(EnsDb.Mmusculus.v79)
edb <- EnsDb.Mmusculus.v79

# add gene symbols 
volcano_prep <- function(res) {
  df <- as.data.frame(res) |>
    tibble::rownames_to_column("ensembl_id")

  df$symbol <- mapIds(x = edb,
                      keys = df$ensembl_id,
                      column = "SYMBOL",
                      keytype = "GENEID")

  df <- df[!is.na(df$symbol), ]
  df <- df[order(df$padj), ]
  df <- df[!duplicated(df$symbol), ]
  return(df)
}

## MALES
volcano_M_WTvsDP16_hyp <- volcano_prep(res_M_WTvsDP16_hyp)
volcano_M_WTvs2N_hyp   <- volcano_prep(res_M_WTvs2N_hyp)
volcano_M_DP16vs2N_hyp <- volcano_prep(res_M_DP16vs2N_hyp)

## FEMALES
volcano_F_WTvsDP16_hyp <- volcano_prep(res_F_WTvsDP16_hyp)
volcano_F_WTvs2N_hyp   <- volcano_prep(res_F_WTvs2N_hyp)
volcano_F_DP16vs2N_hyp <- volcano_prep(res_F_DP16vs2N_hyp)

# build the ranked gene lists 
prep_vector_with_symbols <- function(named_vector) {
  original_ids <- names(named_vector)

  gene_symbols <- mapIds(x = edb,
                         keys = original_ids,
                         column = "SYMBOL",
                         keytype = "GENEID",
                         multiVals = "first")

  data.frame(
    Gene_ID = original_ids,
    Gene_Symbol = as.character(gene_symbols),
    Wald_Statistic = as.numeric(named_vector),
    stringsAsFactors = FALSE
  )
}

output_dir <- "results/hypothalamus"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


## males
write.csv(prep_vector_with_symbols(geneList_M_WTvsDP16_hyp),
          file = file.path(output_dir, "GeneList_Male_WT_vs_Dp16_hyp.csv"), row.names = FALSE)

write.csv(prep_vector_with_symbols(geneList_M_WTvs2N_hyp),
          file = file.path(output_dir, "GeneList_Male_WT_vs_APP2N_hyp.csv"), row.names = FALSE)

write.csv(prep_vector_with_symbols(geneList_M_DP16vs2N_hyp),
          file = file.path(output_dir, "GeneList_Male_Dp16_vs_APP2N_hyp.csv"), row.names = FALSE)

## females
write.csv(prep_vector_with_symbols(geneList_F_WTvsDP16_hyp),
          file = file.path(output_dir, "GeneList_Female_WT_vs_Dp16_hyp.csv"), row.names = FALSE)

write.csv(prep_vector_with_symbols(geneList_F_WTvs2N_hyp),
          file = file.path(output_dir, "GeneList_Female_WT_vs_APP2N_hyp.csv"), row.names = FALSE)

write.csv(prep_vector_with_symbols(geneList_F_DP16vs2N_hyp),
          file = file.path(output_dir, "GeneList_Female_Dp16_vs_APP2N_hyp.csv"), row.names = FALSE)