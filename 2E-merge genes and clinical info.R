# R script: 2E-merge genes and clinical info.R
# Author: Inge Seim
# Use this script to set your parameters

# one can load the gene expression data set from `B-load TOIL TCGA-TARGET-GTex normalised data.R` 
load("TOIL.all.genes.Robj") 

# ######################################################################################
# let us merge the gene expression and clinical data and keep only e.g. PRAD samples
# ######################################################################################
mergedDF <- merge(clinicalDF, geneDF, by="row.names", all=FALSE)  # merge by row names (by=0 or by="row.names")
#@ head(mergedDF) 
#
# colnames(mergedDF) 
# [1] "Row.names"       "Sample.ID"       "Sample.Type"     "DFS"             "DFS.EVENT"            
mergedDF[,c("Row.names","Sample.ID")] # these are all the same. Kill off the first two 
mergedDF <- subset(mergedDF, select=-c(Row.names))
names(mergedDF)[names(mergedDF)=="Sample.ID"] <- "sample"
names(mergedDF)[names(mergedDF)=="DFS.EVENT"] <- "DFS_EVENT"
#
save(mergedDF,file=paste("TCGA-",dataset,".all.genes.and.clinicalRobj",sep=""))
# size=57.9MB
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#
# NOW HAVE A DATA FRAME CONTAINING THE GENE EXPRESSION AND CLINICAL DATA FOR E.G. TCGA-PRAD
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 

# ~~~~~~~~~~~~~~~~~~~~~~~~ THIS IS THE END ~~~~~~~~~~~~~~~~~~~~~~~~ #

