# R script: 2A-set parameters.R
# Author: Inge Seim
# Use this script to set your parameters

# remove all objects
rm(list = ls(all = TRUE))

# load the required libraries
library(data.table) # to load very large data files
library(survival)
library(ggplot2)
library(ggthemes)
library(survminer) # required for ggsurvplot

# go to working directory 
setwd("~/Downloads/R")

# set data set of interest
dataset <- "PRAD"  # TCGA prostate cancer  
# ~~~~~~~~~~~~~~~~~~~~~~~~ THIS IS THE END ~~~~~~~~~~~~~~~~~~~~~~~~ #
