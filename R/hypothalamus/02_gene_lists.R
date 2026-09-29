#results x6 diff results 

##MALE
#WT vs Dp16
res_M_WTvsDP16_hyp <- results(dds_hyp, contrast=list(c("Genotype_Dp16_vs_WT", "SexMale.GenotypeDp16")))
geneList_M_WTvsDP16_hyp <- res_M_WTvsDP16_hyp$stat
names(geneList_M_WTvsDP16_hyp) <- rownames(res_M_WTvsDP16_hyp)
geneList_M_WTvsDP16_hyp <- geneList_M_WTvsDP16_hyp[!is.na(geneList_M_WTvsDP16_hyp)]
geneList_M_WTvsDP16_hyp <- geneList_M_WTvsDP16_hyp[names(geneList_M_WTvsDP16_hyp) != ""]
geneList_M_WTvsDP16_hyp <- sort(geneList_M_WTvsDP16_hyp, decreasing=TRUE)

#WT vs APP2N
res_M_WTvs2N_hyp <- results(dds_hyp, contrast=list(c("Genotype_APP2N_vs_WT", "SexMale.GenotypeAPP2N")))
geneList_M_WTvs2N_hyp <- res_M_WTvs2N_hyp$stat
names(geneList_M_WTvs2N_hyp) <- rownames(res_M_WTvs2N_hyp)
geneList_M_WTvs2N_hyp <- geneList_M_WTvs2N_hyp[!is.na(geneList_M_WTvs2N_hyp)]
geneList_M_WTvs2N_hyp <- geneList_M_WTvs2N_hyp[names(geneList_M_WTvs2N_hyp) != ""]
geneList_M_WTvs2N_hyp <- sort(geneList_M_WTvs2N_hyp, decreasing=TRUE)

#Dp16 vs APP2N
res_M_DP16vs2N_hyp <- results(dds_hyp, contrast=list(c("Genotype_Dp16_vs_WT","SexMale.GenotypeDp16"),
                                             c("Genotype_APP2N_vs_WT","SexMale.GenotypeAPP2N")))
geneList_M_DP16vs2N_hyp <- res_M_DP16vs2N_hyp$stat
names(geneList_M_DP16vs2N_hyp) <- rownames(res_M_DP16vs2N_hyp)
geneList_M_DP16vs2N_hyp <- geneList_M_DP16vs2N_hyp[!is.na(geneList_M_DP16vs2N_hyp)]
geneList_M_DP16vs2N_hyp <- geneList_M_DP16vs2N_hyp[names(geneList_M_DP16vs2N_hyp) != ""]
geneList_M_DP16vs2N_hyp <- sort(geneList_M_DP16vs2N_hyp, decreasing=TRUE)


##FEMALE

#WT vs Dp16 
res_F_WTvsDP16_hyp <- results(dds_hyp, name="Genotype_Dp16_vs_WT")
geneList_F_WTvsDP16_hyp <- res_F_WTvsDP16_hyp$stat
names(geneList_F_WTvsDP16_hyp) <- rownames(res_F_WTvsDP16_hyp)
geneList_F_WTvsDP16_hyp <- sort(geneList_F_WTvsDP16_hyp, decreasing=TRUE)
geneList_F_WTvsDP16_hyp <- geneList_F_WTvsDP16_hyp[!is.na(geneList_F_WTvsDP16_hyp)]
geneList_F_WTvsDP16_hyp <- geneList_F_WTvsDP16_hyp[names(geneList_F_WTvsDP16_hyp) != ""]

#WT vs APP2N
res_F_WTvs2N_hyp <- results(dds_hyp, name="Genotype_APP2N_vs_WT")
geneList_F_WTvs2N_hyp <- res_F_WTvs2N_hyp$stat
names(geneList_F_WTvs2N_hyp) <- rownames(res_F_WTvs2N_hyp)
geneList_F_WTvs2N_hyp <- sort(geneList_F_WTvs2N_hyp, decreasing=TRUE)
geneList_F_WTvs2N_hyp <- geneList_F_WTvs2N_hyp[!is.na(geneList_F_WTvs2N_hyp)]
geneList_F_WTvs2N_hyp <- geneList_F_WTvs2N_hyp[names(geneList_F_WTvs2N_hyp) != ""]

#Dp16 vs APP2N
res_F_DP16vs2N_hyp   <- results(dds_hyp, contrast=list(c("Genotype_Dp16_vs_WT"),
                                               c("Genotype_APP2N_vs_WT")))
geneList_F_DP16vs2N_hyp <- res_F_DP16vs2N_hyp$stat
names(geneList_F_DP16vs2N_hyp) <- rownames(res_F_DP16vs2N_hyp)
geneList_F_DP16vs2N_hyp <- geneList_F_DP16vs2N_hyp[!is.na(geneList_F_DP16vs2N_hyp)]
geneList_F_DP16vs2N_hyp <- geneList_F_DP16vs2N_hyp[names(geneList_F_DP16vs2N_hyp) != ""]
geneList_F_DP16vs2N_hyp <- sort(geneList_F_DP16vs2N_hyp, decreasing=TRUE)