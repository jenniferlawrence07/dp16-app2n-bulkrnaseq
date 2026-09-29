#results x6 diff results 

##MALE
#WT vs Dp16
res_M_WTvsDP16 <- results(dds, contrast=list(c("Genotype_Dp16_vs_WT", "SexMale.GenotypeDp16")))
geneList_M_WTvsDP16 <- res_M_WTvsDP16$stat
names(geneList_M_WTvsDP16) <- rownames(res_M_WTvsDP16)
geneList_M_WTvsDP16 <- geneList_M_WTvsDP16[!is.na(geneList_M_WTvsDP16)]
geneList_M_WTvsDP16 <- geneList_M_WTvsDP16[names(geneList_M_WTvsDP16) != ""]
geneList_M_WTvsDP16 <- sort(geneList_M_WTvsDP16, decreasing=TRUE)

#WT vs APP2N
res_M_WTvs2N <- results(dds, contrast=list(c("Genotype_APP2N_vs_WT", "SexMale.GenotypeAPP2N")))
geneList_M_WTvs2N <- res_M_WTvs2N$stat
names(geneList_M_WTvs2N) <- rownames(res_M_WTvs2N)
geneList_M_WTvs2N <- geneList_M_WTvs2N[!is.na(geneList_M_WTvs2N)]
geneList_M_WTvs2N <- geneList_M_WTvs2N[names(geneList_M_WTvs2N) != ""]
geneList_M_WTvs2N <- sort(geneList_M_WTvs2N, decreasing=TRUE)

#Dp16 vs APP2N
res_M_DP16vs2N <- results(dds, contrast=list(c("Genotype_Dp16_vs_WT","SexMale.GenotypeDp16"),
                                             c("Genotype_APP2N_vs_WT","SexMale.GenotypeAPP2N")))
geneList_M_DP16vs2N <- res_M_DP16vs2N$stat
names(geneList_M_DP16vs2N) <- rownames(res_M_DP16vs2N)
geneList_M_DP16vs2N <- geneList_M_DP16vs2N[!is.na(geneList_M_DP16vs2N)]
geneList_M_DP16vs2N <- geneList_M_DP16vs2N[names(geneList_M_DP16vs2N) != ""]
geneList_M_DP16vs2N <- sort(geneList_M_DP16vs2N, decreasing=TRUE)


##FEMALE

#WT vs Dp16 
res_F_WTvsDP16 <- results(dds, name="Genotype_Dp16_vs_WT")
geneList_F_WTvsDP16 <- res_F_WTvsDP16$stat
names(geneList_F_WTvsDP16) <- rownames(res_F_WTvsDP16)
geneList_F_WTvsDP16 <- sort(geneList_F_WTvsDP16, decreasing=TRUE)
geneList_F_WTvsDP16 <- geneList_F_WTvsDP16[!is.na(geneList_F_WTvsDP16)]
geneList_F_WTvsDP16 <- geneList_F_WTvsDP16[names(geneList_F_WTvsDP16) != ""]

#WT vs APP2N
res_F_WTvs2N <- results(dds, name="Genotype_APP2N_vs_WT")
geneList_F_WTvs2N <- res_F_WTvs2N$stat
names(geneList_F_WTvs2N) <- rownames(res_F_WTvs2N)
geneList_F_WTvs2N <- sort(geneList_F_WTvs2N, decreasing=TRUE)
geneList_F_WTvs2N <- geneList_F_WTvs2N[!is.na(geneList_F_WTvs2N)]
geneList_F_WTvs2N <- geneList_F_WTvs2N[names(geneList_F_WTvs2N) != ""]

#Dp16 vs APP2N
res_F_DP16vs2N   <- results(dds, contrast=list(c("Genotype_Dp16_vs_WT"),
                                               c("Genotype_APP2N_vs_WT")))
geneList_F_DP16vs2N <- res_F_DP16vs2N$stat
names(geneList_F_DP16vs2N) <- rownames(res_F_DP16vs2N)
geneList_F_DP16vs2N <- geneList_F_DP16vs2N[!is.na(geneList_F_DP16vs2N)]
geneList_F_DP16vs2N <- geneList_F_DP16vs2N[names(geneList_F_DP16vs2N) != ""]
geneList_F_DP16vs2N <- sort(geneList_F_DP16vs2N, decreasing=TRUE)
