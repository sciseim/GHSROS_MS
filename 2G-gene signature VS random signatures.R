# R script: 2G-gene signature VS random signatures.R
# Author: Inge Seim
# Use this script to compare your gene expression signature to N random signatures
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# load gene expression data and subset to fit your data set
# here: TCGA-PRAD
all.TCGA.genes <- tempdata[,5:length(colnames(tempdata))] # 489 PRAD-TCGA SAMPLES AND 58,581 genes
row.names(all.TCGA.genes) <-tempdata$sample
# length(all.TCGA.genes) # 57581 genes
# colnames(all.TCGA.genes) # genes only


# what is the length or your signature?
length(genesignature) # e.g. 34-gene signature
dataset <- "PRAD" # e.g. TCGA-PRAD


# prepare for k-means clustering and survival analysis of N random genes of the same size as your gene signature
inputDF <- all.TCGA.genes 
save(inputDF,file=paste("RandomGeneSignatureCheckInput-TCGA-",dataset,"-",length(genesignature),"-gene signature.Robj",sep="")
) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set parameters
signaturelength <- length(genesignature) # e.g. 34
ITERATIONS <- 100000 # should be minumum 1000
#
# k-means settings
set.seed(20)
k <- 2 # we want *k* groups to be 2. i.e. High and Low expression of the genes in the signature
numberofruns <- 500
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#  RUN
system.time(source("random_signature_kmeans_survival_loop.R"))


# ##############################################################
# SAVE THE OUTPUT
# ##############################################################
signature.analysis.DF <- do.call(rbind.data.frame, signature.analysis.list)  # convert the list to a data frame!
#
# save the output data
save(signature.analysis.DF,file=paste("signature.analysis-",ITERATIONS,"-random.genes.Robj",sep=""))


# ################################################################################################################################
# draw a density plot 
# based style from 'SigCheck' and Extended Data Figure 5 of Nature.2015 16;523(7560):313-7 ('SigCheck' author co-author)
# 'SigCheck' does not implement k-means, however.
genesignatureDFS.p.value <- signaturePvalue  # the DFS P of our 34-gene signature 
# genesignatureDFS.p.value <- 0.00019 # manual
randomscores <- as.matrix(signature.analysis.DF) # the DFS P of random signatures of the same size
#
class(randomscores) # matrix
NumberOfRandomSignaturesTested <- length(randomscores) # e.g. 1000  
NumberOfRandomSignaturesSmallerOrEqualToOurSignature <- sum(randomscores <= genesignatureDFS.p.value) 
empirical.P.value <- (NumberOfRandomSignaturesSmallerOrEqualToOurSignature / NumberOfRandomSignaturesTested)
Percentile <- 1-empirical.P.value

titlestring <- "Survival: Random Signature"
subtitle <- paste(Percentile," (Tests: ",NumberOfRandomSignaturesTested," P=",empirical.P.value,")",sep="")                           
LOGrandomscores <- log10(randomscores)
LOG005 <- log10(0.05) # -1.30103
LOGgenesig <- log10(genesignatureDFS.p.value) # e.g. -3.720787


pdf(paste(ITERATIONS,"-RANDOM-gene-signatures-density-plot.pdf",sep=""))
p <- plot(density(LOGrandomscores),main=titlestring,sub=subtitle)
abline(v=LOGgenesig,col="#1c6cab",lwd=3)
abline(v=LOG005,col="#1c6cab",lty="dotted",lwd=2)
legend("topright",legend=c(sprintf("Signature  P-val %0.5f",
                                   genesignatureDFS.p.value),sprintf("Significant P-val %0.2f",
                                                                     0.05)),col="#1c6cab",lty=c("solid","dotted"),lwd=c(3,2))
print(p)
dev.off()
# ~~~~~~~~~~~~~~~~~~~~~~~~ THIS IS THE END ~~~~~~~~~~~~~~~~~~~~~~~~ #

