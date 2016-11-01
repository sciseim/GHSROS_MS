# R script: 1A-SCAN.UPC.R
# Author: Inge Seim
# This script will run the Single-channel array normalization (SCAN) and Universal exPression Codes (UPC) methods on Affymetrix CEL files

# remove all objects
rm(list = ls(all = TRUE))

# load the required libraries
library(SCAN.UPC)
library(pd.huex.1.0.st.v2)

# go to working directory with data sets in different directories
# e.g. ./GSE000001 has 10 CEL files; ./GSE999999 has 5 CEL files
setwd("/Volumes/ExonArray/GHSROS/data/CEL")
startdirectory <- getwd()

# list directories
CELdirs <- list.dirs('.', recursive=FALSE,full.names=FALSE)
length(CELdirs)

# loop
for(i in 1:length(CELdirs))
{   # START LOOP
dataset <- CELdirs[i]
setwd(paste(startdirectory,"/",CELdirs[i],sep=""))

# run SCAN
normalisedname <- paste("normalised.",dataset,sep="")
normalisedname <- SCAN("*.CEL", outFilePath = NA, convThreshold = 0.01, annotationPackageName = "pd.huex.1.0.st.v2", probeSummaryPackage = NA, probeLevelOutDirPath = NA, exonArrayTarget="probeset", batchFilePath=NA, verbose = TRUE)
# annotationPackageName = "pd.huex.1.0.st.v2"
# default is 'annotationPackageName = NA'. Add 'pd.huex.1.0.st.v2' to force new probe set annotation for older Exon Arrays.

# run UPC
UPCnormalisedname <- paste("UPCnormalised.",dataset,sep="")
UPCnormalisedname <- UPC_Generic_ExpressionSet(normalisedname, sequenceFeatureName = NA, modelType = "nn", convThreshold = 0.001, higherValuesIndicateHigherExpression = TRUE, verbose = TRUE)

# convert the ExpressionSet objects to data frames
normalisednameDF <- paste("normalisedDF.",dataset,sep="")
UPCnormalisednameDF <- paste("UPCnormalisedDF.",dataset,sep="")
normalisednameDF  <- as.data.frame(exprs(normalisedname))
UPCnormalisednameDF <- as.data.frame(exprs(UPCnormalisedname))

# save the data (large files, so compress using RDS)
# SCAN
RDSnameN <- paste("normalised.",dataset,".rds",sep="")
RDSnameNdf <- paste("normalisedDF.",dataset,".rds",sep="")
saveRDS(normalisedname,RDSnameN) # SCAN ExpressionSet Object
saveRDS(normalisednameDF,file=RDSnameNdf) # SCAN data frame
# UPC
RDSnameUPC <- paste("UPC.",dataset,".rds",sep="")
RDSnameUPCdf <- paste("UPCDF.",dataset,".rds",sep="")
saveRDS(UPCnormalisedname,file=RDSnameUPC) # UPC ExpressionSet Object
saveRDS(UPCnormalisednameDF,file=RDSnameUPCdf) # UPC data frame

} # END LOOP
# ~~~~~~~~~~~~~~~~~~~~~~~~ THIS IS THE END ~~~~~~~~~~~~~~~~~~~~~~~~ #
