# R script: 2C-input gene signature.R
# Author: Inge Seim
# Use this script to input your gene expression signature of interest

# 34-gene singature from our manuscript
genesignature <- c("FBXL16","DIRAS1","CRIP2","TP53I11","MUC5B","TFF2","ZNF467","NUDT11","PARM1","CNTN1","ST6GAL1","UNC80","EYA1","MUM1L1","HSPB8","VWA5A","KLF9","DMD","CHRDL1","MCTP2","RYR2","ANGPT1","NUDT10","STOX2","CAPN6","EGF","RNASEL","AASS","HEPH","LRCH2","LCP1","CPA6","IFI16","LRIG1")
class(genesignature) # character
# or, if one prefers to load a text file
genesignature <- readLines("./data/the_signature.txt")
class(genesignature) # character


# Save a data.frame that only contain the gene expression data for your gene signature
# Useful for e.g. drawing heat maps
genesignatureDF <- geneDF
colnames(genesignatureDF) # genes
genesignatureDF <- genesignatureDF[,genesignature]  # only keep the gene signature genes
#@ genesignatureDF["TCGA-2A-A8VV-01","FBXL16"] # checked vs xena online: CORRECT
saveRDS(genesignatureDF,file="genesignatureDF.Rds")
# ~~~~~~~~~~~~~~~~~~~~~~~~ THIS IS THE END ~~~~~~~~~~~~~~~~~~~~~~~~ #
