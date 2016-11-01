# R script: 2F-gene signature k-means clustering and Kaplan-Meier survival analysis.R
# Author: Inge Seim
# Perform k-means clustering to stratify your patients, followed by KM survival analysis


colnames(mergedDF) # genes
# get rid of NA DFS values. Will 'choke' some R packages (e.g. 'SigCheck')
# need a format compatible with the 'survival' R package
tempdata <- mergedDF

# omit rows with DFS or OS <NA>
# here, we are working with TCGA-PRAD (no OS values)
# DFS
# class(tempdata$DFS_EVENT) # factor
tempdata$DFS_EVENT <- as.character(tempdata$DFS_EVENT)
tempdata <- tempdata[!is.na(tempdata$DFS_EVENT),]
tempdata[tempdata=="DiseaseFree"]<-0
tempdata[tempdata=="Recurred/Progressed"]<-1
tempdata$DFS_EVENT <- as.numeric(tempdata$DFS_EVENT)
#
# OS
# class(tempdata$OS_EVENT) # factor
# tempdata$OS_EVENT <- as.character(tempdata$OS_EVENT)
# tempdata <- tempdata[!is.na(tempdata$OS_EVENT),]
# tempdata[tempdata=="LIVING"]<-0
# tempdata[tempdata=="DECEASED"]<-1
# tempdata$OS_EVENT <- as.numeric(tempdata$OS_EVENT)

tempdataGenesOnly <- tempdata[,genesignature] # keep your gene expression signature of interest
# *ONLY* GENE EXPRESSION DATA HERE
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 


# ##############################################################
# k-means clustering
# ##############################################################
# sort by gene name
# class(names(testDF))
#@ colnames(tempdataGenesOnly) # only gene expression data here!

set.seed(20)
k <- 2 # we want *k* groups to be 2. i.e. High and Low expression of the genes in the signature
numberofruns <- 500
samplecluster <- kmeans(tempdataGenesOnly, k, nstart = numberofruns)


# assign
samplecluster$cluster <- as.factor(samplecluster$cluster)

# ADD THE GROUP INFORMATION TO YOUR ORIGINAL DATA FRAME
# append cluster assignment
tempdataGenesOnly <- data.frame(tempdataGenesOnly, samplecluster$cluster)
# rename to fit preferred style
names(tempdataGenesOnly)[names(tempdataGenesOnly) == 'samplecluster.cluster'] <- 'group'

# assign a more useful name
levels(tempdataGenesOnly$group)[levels(tempdataGenesOnly$group)=="1"] <- "cluster1"
levels(tempdataGenesOnly$group)[levels(tempdataGenesOnly$group)=="2"] <- "cluster2"


# add DFS and DFS_EVENT data from tempdata above
tempdataGenesOnly$DFS_EVENT <- as.numeric(tempdata$DFS_EVENT)
tempdataGenesOnly$DFS <- (tempdata$DFS)


# ##############################################################
# survival analysis
# ##############################################################
survdiff(Surv(DFS, DFS_EVENT) ~ group ,data=tempdataGenesOnly , rho=0)  

DFSsurv <- survdiff(Surv(DFS, DFS_EVENT) ~ group ,data=tempdataGenesOnly, rho=0)
p.val.DFS <- 1 - pchisq(DFSsurv$chisq, length(DFSsurv$n) - 1)
signaturePvalue <- p.val.DFS

# save output
tempdataYear <- tempdataGenesOnly
tempdataYear$DFS # months
tempdataYear$OS # months
tempdataYear$DFS <- tempdataYear$DFS/12 # years
tempdataYear$OS <- tempdataYear$OS/12 # years
pdf(paste(dataset,"-kmeans-survival-DFS-UCSC-TOIL-normalised-gene-counts.pdf",sep=""))
fit=survfit(Surv(DFS, DFS_EVENT) ~ group,data=tempdataYear)
p <- ggsurvplot(fit, risk.table = FALSE,
                pval = TRUE, risk.table.y.text.col = FALSE, break.time.by=1,conf.int=FALSE, censor=TRUE, legend="top", palette=c("#f42732","#000000"))
print(p)
dev.off()
# ~~~~~~~~~~~~~~~~~~~~~~~~ THIS IS THE END ~~~~~~~~~~~~~~~~~~~~~~~~ #
