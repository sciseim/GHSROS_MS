# R script: 2H-gene signature VS known signatures.R
# Author: Inge Seim
# Use this script to compare your gene expression signature and known signatures
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
PCaSigTitles <- c("Bibikova-2007","Peng-2014","Glinsky-2005",	"Thompson-2012", "Li-2013",	"Penney-2011",	"Rajan-2014",	"Chen-2011",	"Bismar-2006",	"Nakagawa-2008", "Singh-2002",	"Lapointe-2004",	"Ramaswamy-2003",	"Glinsky-2004",	"Wu-2013",	"Chen-2012",	"Zhao-2016",	"Yu-2007",	"Prensner-SChLAP1",	"Varambally-2005",	"Saal-2007",	"Cuzick-2011")
PCaSigLogRankP <- c("0.000155","0.0367","1.05E-05","0.217","0.00192","0.000737","0.03","0.144","0.557","0.00179","0.332","0.00105","0.0161","0.0783","0.000172","0.222","0.042","0.193","0.0757","0.968","1.37E-07","1.53E-07")
PCaSignaturesDF <- data.frame(PCaSigTitles, PCaSigLogRankP)


# ################################################################################################################################
# draw a density plot 
# based style from 'SigCheck' and Extended Data Figure 5 of Nature.2015 16;523(7560):313-7 ('SigCheck' author co-author)
# 'SigCheck' does not implement k-means, however.

genesignatureDFS.p.value <- signaturePvalue  # the DFS P of our 34-gene signature 
# genesignatureDFS.p.value <- 0.00019 # manual

PCaSigLogRankP <- as.numeric(PCaSigLogRankP) # 

NumberOfSignaturesTested <- length(PCaSigLogRankP) # e.g. 22
NumberOfSignaturesSmallerOrEqualToOurSignature <- sum(PCaSigLogRankP <= genesignatureDFS.p.value)  # e.g. 5
empirical.P.value <- (NumberOfSignaturesSmallerOrEqualToOurSignature / NumberOfSignaturesTested)
Percentile <- 1-empirical.P.value
#
titlestring <- "Survival: Known Signatures"
subtitle <- paste(Percentile," (Tests: ",NumberOfSignaturesTested," p=",empirical.P.value,")",sep="")                           


LOGPCaSigLogRankP <- log10(PCaSigLogRankP)
LOG005 <- log10(0.05) # -1.30103
LOGgenesig <- log10(genesignatureDFS.p.value) # 

pdf(paste("KNOWN-gene-signatures-density-plot.pdf",sep=""))
p <- plot(density(LOGPCaSigLogRankP),main=titlestring,sub=subtitle)
abline(v=LOGgenesig,col="#1c6cab",lwd=3)
abline(v=LOG005,col="#1c6cab",lty="dotted",lwd=2)
legend("topright",legend=c(sprintf("Signature  P-val %0.5f",
                                   genesignatureDFS.p.value),sprintf("Significant P-val %0.2f",
                                                                     0.05)),col="#1c6cab",lty=c("solid","dotted"),lwd=c(3,2))
print(p)
dev.off()
# ~~~~~~~~~~~~~~~~~~~~~~~~ THIS IS THE END ~~~~~~~~~~~~~~~~~~~~~~~~ #
