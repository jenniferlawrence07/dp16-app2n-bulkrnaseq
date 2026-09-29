
install.packages('RMySQL', repos='http://cran.us.r-project.org')
install.packages('BiocManager', repos='http://cran.us.r-project.org')
install.packages('Rsubread', repos='http://cran.us.r-project.org')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")


files=c("Hoeffer_101_S16_L005_Aligned.out.bam",
        "Hoeffer_118_S17_L005_Aligned.out.bam",
        "Hoeffer_123_S18_L005_Aligned.out.bam",
        "Hoeffer_125_S32_L005_Aligned.out.bam",
        "Hoeffer_126_S36_L005_Aligned.out.bam",
        "Hoeffer_28_S31_L005_Aligned.out.bam",
        "Hoeffer_31_S19_L005_Aligned.out.bam",
        "Hoeffer_33_S25_L005_Aligned.out.bam",
        "Hoeffer_37_S13_L005_Aligned.out.bam",
        "Hoeffer_39_S26_L005_Aligned.out.bam",
        "Hoeffer_40_S27_L005_Aligned.out.bam",
        "Hoeffer_41_S29_L005_Aligned.out.bam",
        "Hoeffer_44_S24_L005_Aligned.out.bam",
        "Hoeffer_47_S33_L005_Aligned.out.bam",
        "Hoeffer_48_S34_L005_Aligned.out.bam",
        "Hoeffer_52_S21_L005_Aligned.out.bam",
        "Hoeffer_54_S22_L005_Aligned.out.bam",
        "Hoeffer_55_S14_L005_Aligned.out.bam",
        "Hoeffer_56_S23_L005_Aligned.out.bam",
        "Hoeffer_74_S30_L005_Aligned.out.bam",
        "Hoeffer_75_S28_L005_Aligned.out.bam",
        "Hoeffer_77_S15_L005_Aligned.out.bam",
        "Hoeffer_82_S35_L005_Aligned.out.bam",
        "Hoeffer_88_S20_L005_Aligned.out.bam")

Gtf=("/pl/archive/Link_Lab_archive/griffioen/Mus_musculus.GRCm39.111.gtf")

###Unique
uniqueP=featureCounts(files,
isGTFAnnotationFile = TRUE,
annot.ext = gtf,
GTF.attrType = "gene_id",
GTF.featureType = "gene",
useMetaFeatures = FALSE,
allowMultiOverlap = FALSE,
countMultiMappingReads = FALSE,
isPairedEnd=TRUE,
nthreads = 8,
strandSpecific = 1)

write.table(x=data.frame(
uniqueP$annotation[,c("GeneID")],
uniqueP$counts,stringsAsFactors=FALSE),
file="APPDp16Hyp_counts.txt",
quote=FALSE,sep="\t",row.names=FALSE)
