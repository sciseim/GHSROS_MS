# R script: 2D-clinical information.R
# Author: Inge Seim
# This script will load TCGA clinical information from cBioPortal (continuously-updated)

# ######################################################################################
# load clinical information
# ######################################################################################
# download from the cBioPortal for Cancer Genomics (continously updated)
#   J. Gao et al., Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. Sci Signal 6, pl1 (2013).
#.	E. Cerami et al., The cBio cancer genomics portal: an open platform for exploring multidimensional cancer genomics data. Cancer Discov 2,  #   401-404 (2012).
#
# here, we also load a list of the samples used to assess
system.time(clinicalDF <- fread(paste("./data/cBioPortal/",dataset,".txt",sep=""), header=T, sep="\t", stringsAsFactors=T,showProgress=TRUE)) # RNA
clinicalDF <- as.data.frame(clinicalDF)
names(clinicalDF) <- gsub(" ", ".", names(clinicalDF)) # space incompatible
names(clinicalDF) <- gsub("\\(", ".", names(clinicalDF)) # space incompatible
names(clinicalDF) <- gsub("\\)", ".", names(clinicalDF)) # space incompatible

row.names(clinicalDF) <- clinicalDF$Sample.ID

# rename DFS events to something that is easier to read
# [53] "Overall.Survival..Months."                                                                  
# [54] "Overall.Survival.Status"  
# [24] "Disease.Free..Months."                                                                      
# [25] "Disease.Free.Status" 
names(clinicalDF)[names(clinicalDF) == 'Disease.Free..Months.'] <- 'DFS'
names(clinicalDF)[names(clinicalDF) == 'Disease.Free.Status'] <- 'DFS.EVENT'
clinicalDF <- clinicalDF[,c("Sample.ID","Sample.Type","DFS","DFS.EVENT")] # clean up

# ~~~~~~~~~~~~~~~~~~~~~~~~ THIS IS THE END ~~~~~~~~~~~~~~~~~~~~~~~~ #
