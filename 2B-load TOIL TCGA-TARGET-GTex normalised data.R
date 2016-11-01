# R script: 2B-load TOIL TCGA-TARGET-GTex normalised data.R
# Author: Inge Seim
# This script will load normalised gene expression data from ~60,000 genes and 20,000 samples

# ######################################################################################
# load gene expression
# ######################################################################################
#
# file with normalised gene expression values. Here, normalised gene counts from the UCSC TOIL RNA-seq recompute
# J. Vivian et al., Rapid and efficient analysis of 20,000 RNA-seq samples with Toil. bioRxiv,  (2016).
#
# 1. download log2(norm_counts+1) from TCGA TARGET GTEx gene expression by UCSC TOIL RNA-seq recompute (see https://goo.gl/u8acdp)
#       wget https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count;
#    Reprocessed GTEx, TCGA and TARGET with the same pipeline, removing computational batch effects 
#    GTEx, TCGA and TARGET all use different expression calling programs. We used the same program for all of them.
#    log2(norm_counts+1), mean is subtracted per column across 19,120 samples. 
# 2. ~8GB file, so transpose in the UNIX command prompt using Jim Kent's rowsToCols [https://github.com/ENCODE-DCC/kentUtils]  
#       rowsToCols TcgaTargetGtex_RSEM_Hugo_norm_count TcgaTargetGtex_RSEM_Hugo_norm_count_ROWS2COL.txt;
# 
system.time(geneDF <- fread("TcgaTargetGtex_RSEM_Hugo_norm_count_ROWS2COL.txt", header=T, sep="\t", stringsAsFactors=F,showProgress=TRUE)) # RNA-SEQ GENE EXP
# size = 8.02 GB
geneDF <- as.data.frame(geneDF)

# # These are HUGO gene names, as the file name suggests. 
# Annotation: Gencode V23 comprehensive annotation (CHR) http://www.gencodegenes.org/releases/23.html
# row.names(geneDF) # empty
geneDF <- subset(geneDF, !duplicated(geneDF[,1]))  # the current version has a couple of duplicated sample names
# ~20k samples. Will take some time to run.

geneDF$sample # samples
row.names(geneDF) <- geneDF$sample
# kill sample column
geneDF <- geneDF[,-1]  # subset(geneDF, select=-c(sample)) too slow. Will search 20k names
colnames(geneDF) # genes
row.names(geneDF # samples
head(geneDF[1:2,1:5]) # all 19,120 samples
#                     sample             CTD-2588J6.1 RP11-433M22.1
# TCGA.OR.A5JX.01 TCGA-OR-A5JX-01            0             0
# TCGA.HV.A5A5.01 TCGA-HV-A5A5-01            0             0

# save Robject
save(geneDF,file="TOIL.all.genes.Robj")
#@ load("TOIL.all.genes.Robj") # to load later
# size = 1.96GB 
# ~~~~~~~~~~~~~~~~~~~~~~~~ THIS IS THE END ~~~~~~~~~~~~~~~~~~~~~~~~ #
