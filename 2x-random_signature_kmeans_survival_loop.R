# R script: 2x-random_signature_kmeans_survival_loop.R
# Author: Inge Seim

remove(signature.analysis.list)
signature.analysis.list <- list()


for (i in 1:ITERATIONS)
{   # START LOOP
  temp <- sample(inputDF,size=signaturelength, replace=FALSE)
  
  #@ print(colnames(temp))
  #@ length(temp) # e.g. 34
  
  # temp2 <- temp
  samplecluster <- kmeans(temp, k, nstart = numberofruns)
  # assign
  samplecluster$cluster <- as.factor(samplecluster$cluster)
  # ADD THE GROUP INFORMATION TO YOUR ORIGINAL DATA FRAME
  # append cluster assignment
  temp <- data.frame(temp, samplecluster$cluster)
  # rename to fit preferred style
  names(temp)[names(temp) == 'samplecluster.cluster'] <- 'group'
  # assign a more useful name
  levels(temp$group)[levels(temp$group)=="1"] <- "cluster1"
  levels(temp$group)[levels(temp$group)=="2"] <- "cluster2"
  # add DFS and DFS.EVENT
  temp$DFS_EVENT <- as.character(tempdata$DFS_EVENT)
  temp$DFS <- as.vector(tempdata$DFS)
  # prep
  temp[temp=="DiseaseFree"]<-0
  temp[temp=="Recurred/Progressed"]<-1
  temp$DFS <- as.numeric(temp$DFS )
  temp$DFS_EVENT <- as.numeric(temp$DFS_EVENT )
  # temp$DFS_EVENT <- as.character(temp$DFS_EVENT)
  temp$DFS <- as.vector(temp$DFS)
  
  # ##############################################################
  # survival analysis
  # ##############################################################
  survdiff(Surv(DFS, DFS_EVENT) ~ group ,data=temp , rho=0)  
  # p= 0.00019 
  # yes, the correct data pulled out of the ExpressionSet. ExpressionSet OK to use for SigCheck.
  DFSsurv <- survdiff(Surv(DFS, DFS_EVENT) ~ group ,data=temp, rho=0)
  p.val.DFS <- 1 - pchisq(DFSsurv$chisq, length(DFSsurv$n) - 1)
  
  signature.analysis.list[[i]] <- p.val.DFS  # save for each i (here: TCGA data set)
  
  
  
  # save data for all
#@  print(colnames(temp))
  print(i)
#@  print(signature.analysis.list[[i]])
  # remove(temp)
} # END LOOP
