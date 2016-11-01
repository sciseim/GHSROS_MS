# R script: 4-TCGA-PRAD heatmap.R
# Author: Inge Seim
# Use this script to generate a heatmap of your gene set of interest

# ######################################################################################
# PRAD heat map figure for paper
# ######################################################################################
library(data.table)
library(gplots)
library(RColorBrewer)
library(devtools)

# ######################################################################################
# decide on what data set you want to interrogate
# ######################################################################################
dataset <- "PRAD" # e.g. TCGA-PRAD=prostate cancer



# ######################################################################################
# load gene expression
# ######################################################################################
# note: the gene expression values here were obtained from xena 
# Gene Hugo RSEM norm_counts
# log2(norm_counts+1)

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# LOAD THE DATA FRAME *AFTER* RUNNING `3-survival_analysis_TCGA-datasets.R`
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
load(file=paste("./output/Robjects/",dataset,"ready.Robj",sep="")) # load this since we want k-means cluster group info. too
geneDF <- tempdata
# columns:
# sample (e.g. TCGA-2A-A8VO-01)
# sample_type (e.g. Primary Tumor)
# DFS (e.g. 55.88 months)
# DFS_EVENT (e.g. 0 disease-free or cencored)
# ... gene expression log2(RSEM+1)
# group: what k-means cluster does sample belong to?
#
# if you do not care about cluster etc, load a data frame with sample names and normalised
# gene expressions.

# if you want to add a side bar showing gene expression in another data set, edit ./data/signature-direction.xlsx
# or load a data set with the following structure
# gene	Metastatic
# FBXL16	up
# ....
# NUDT11	down


# ######################################################################################
# load clinical data
# ######################################################################################

system.time(clinicalDF <- fread(paste("./data/cBioPortal/",dataset,".txt",sep=""), header=T, sep="\t", stringsAsFactors=T)) # RNA
clinicalDF <- as.data.frame(clinicalDF) 
head(clinicalDF)
colnames(clinicalDF)
names(clinicalDF) <- gsub(" ", ".", names(clinicalDF)) # space incompatible
names(clinicalDF) <- gsub("\\(", ".", names(clinicalDF)) # space incompatible
names(clinicalDF) <- gsub("\\)", ".", names(clinicalDF)) # space incompatible

# only select a few, common data columns
clinicalDF$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer
clinicalDF$Psa.most.recent.results
clinicalDF$Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Primary.Tumor.T.Stage
# clinicalDF <- subset(clinicalDF, select=c(Sample.ID,Sample.Type,Disease.Free..Months.,Disease.Free.Status,Overall.Survival..Months.,Overall.Survival.Status))
clinicalDF <- subset(clinicalDF, select=c(Sample.ID,Sample.Type,Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer,Psa.most.recent.results,Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Primary.Tumor.T.Stage))

names(clinicalDF)[names(clinicalDF) == 'Sample.ID'] <- 'sample'
names(clinicalDF)[names(clinicalDF) == 'Sample.Type'] <- 'sample_type'
names(clinicalDF)[names(clinicalDF) == 'Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer'] <- 'Gleason.Score'
names(clinicalDF)[names(clinicalDF) == 'Psa.most.recent.results'] <- 'PSA'
names(clinicalDF)[names(clinicalDF) == 'Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Primary.Tumor.T.Stage'] <- 'stage'


# ##############################################################################################################

# ######################################################################################
# merge
# ######################################################################################
mergedDF <- merge(clinicalDF,geneDF, by.x='sample', by.y='sample',all=FALSE) # all=FALSE ... do not keep 
head(mergedDF[1,])
colnames(mergedDF)

# gene expression here
colnames(mergedDF)
geneEXP <- mergedDF
rownames(geneEXP) <- mergedDF$sample
geneEXP <- geneEXP[,11:length(geneEXP)]

head(geneEXP)
geneEXP <- subset(geneEXP, select=-c(group))
# want
# row.name = gene A --sample n
# columns = sample 1 --sample n
# geneEXP <- t(geneEXP) # if you want to have samples as columns
head(geneEXP)
# now have gene names as rows and samples as columns

# working 
value.top <- geneEXP 
# get rid of all NA
ind <- apply(value.top, 1, function(x) all(is.na(x)))
value.top <- value.top[ !ind, ]





# ######################################################################################
# HEAT MAPS
# ######################################################################################

#Load latest version of heatmap.3 function
# source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
source("heatmap.3.R")
#
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

# heatmap colour parameters
cols <- colorRampPalette(c("green", "black", "red"))(n = 256) # classical

# assess the heat map
h <- heatmap.3(scale(value.top), col=cols, scale="column", trace="none",density.info="none", dendrogram="none",Colv=FALSE, key=TRUE,cexRow=0.5, hclustfun=myclust, distfun=mydist, breaks=seq(-3,3,6/256) )
print(h)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ADD SIDE-BAR AND TOP BARS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPACER BAR SPACER BAR SPACER BAR SPACER BAR SPACER BAR SPACER BAR SPACER BAR SPACER BAR
spacercolours=sample(c("white"), length(mergedDF$group), replace = TRUE, prob = NULL)
# SPACER BAR SPACER BAR SPACER BAR SPACER BAR SPACER BAR SPACER BAR SPACER BAR SPACER BAR

# k-MEANS CLUSTER k-MEANS CLUSTER k-MEANS CLUSTER k-MEANS CLUSTER k-MEANS CLUSTER k-MEANS CLUSTER 
groupcolours <- mergedDF$group
class(groupcolours) # character
levels(groupcolours)[levels(groupcolours)=="cluster1"] <- "#f42732"
levels(groupcolours)[levels(groupcolours)=="cluster2"] <- "#000000"
groupcolours <- as.character(groupcolours)
rlab=cbind(groupcolours) # can add as many as we want. E.g. colour for stage
colnames(rlab)=c("group")
# to ensure you are not stuffing up... check out the different cluster sizes
length(which(mergedDF$group == "cluster1")) # 157 samples
length(which(mergedDF$group == "cluster2")) # 332 samples
length(which(mergedDF$group != "cluster1")) # 332 samples
# k-MEANS CLUSTER k-MEANS CLUSTER k-MEANS CLUSTER k-MEANS CLUSTER k-MEANS CLUSTER k-MEANS CLUSTER 

# GLEASON SCORE GLEASON SCORE GLEASON SCORE GLEASON SCORE GLEASON SCORE GLEASON SCORE 
# all Gleasons above a certain score
temp <- mergedDF
temp$Gleason.Score = ifelse(temp$Gleason.Score >= 8, "purple", "grey") # two colours here since two groups
head(temp$Gleason.Score)
Gleasoncolours <- temp$Gleason.Score
rlab=cbind(Gleasoncolours) # can add as many as we want. e.g. colour for stage
colnames(rlab)=c("Gleason")
class(rlab) # matrix
# GLEASON SCORE GLEASON SCORE GLEASON SCORE GLEASON SCORE GLEASON SCORE GLEASON SCORE 


# GRASSO/TAYLOR BAR GRASSO/TAYLOR BAR GRASSO/TAYLOR BAR GRASSO/TAYLOR BAR GRASSO/TAYLOR BAR
# would be nice to indicate whether UP or DOWN
# in Grasso and Taylor Metastasis data sets
library("XLConnect")
excel.file <- file.path("./data/signature-direction.xlsx")
signaturedirection <- readWorksheetFromFile(excel.file, sheet=1)
head(signaturedirection)
signaturedirection <- signaturedirection$Metastatic
class(signaturedirection) # character
# levels(signaturedirection)[levels(signaturedirection)=="up"] <- "red"
# levels(signaturedirection)[levels(signaturedirection)=="down"] <- "green"
head(signaturedirection)

# GRASSO/TAYLOR BAR GRASSO/TAYLOR BAR GRASSO/TAYLOR BAR GRASSO/TAYLOR BAR GRASSO/TAYLOR BAR
column_annotation <- signaturedirection
column_annotation <- gsub("up", "#e80a76", column_annotation)
column_annotation <- gsub("down", "#21823b", column_annotation)
column_annotation <- as.matrix(column_annotation)
colnames(column_annotation) <- c("Metastases Expression")

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# combine
rlab=t(cbind(groupcolours,spacercolours, Gleasoncolours))

#@ pdf("heatmap3-SIDEBARS-scaled-column.pdf")
pdf("TCGA-PRAD-34-gene-signature-heatmap.pdf")
# same as above ... scale by columns = the expression of each gene
h <- heatmap.3(scale(value.top), col=cols, scale="column", trace="none",density.info="none", dendrogram="none",Colv=FALSE, key=TRUE,cexRow=0.5, hclustfun=myclust, distfun=mydist, ColSideColors=column_annotation,RowSideColors=rlab,breaks=seq(-3,3,6/256) )
print(h)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~END~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
