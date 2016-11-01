# R script: 1C-UPC_scatter_plot.R
# Author: Inge Seim
# Description: This script will draw a scatter plot of UPC scores in your data set. Requires output from UPC_parse.R

# remove all objects
rm(list = ls(all = TRUE))

# go to directory
setwd("/Volumes/ExonArray/GHSROS/data/")

# load
mergedDF <- readRDS(file="UPC.mergedDF.rds")

# set a UPC threshold (for drawing a dotted line, does not exclude data points)
UPCthreshold <- 0.10

# decide on probeset of interest and subset
probestofinterestDF <- mergedDF[c("2652604"),] # GHSROS
probestofinterestDF <- as.data.frame(t(probestofinterestDF))  # transpose
names(probestofinterestDF) <- c("probeset") # rename column header to a more generic name

# split into two UPC groups by UPCthreshold
probestofinterestDF$UPCgroup <- as.numeric(probestofinterestDF$probeset) # copy the column
probestofinterestDF$UPCgroup[probestofinterestDF$UPCgroup >= UPCthreshold] <- 1; probestofinterestDF$UPCgroup
probestofinterestDF$UPCgroup[probestofinterestDF$UPCgroup < UPCthreshold] <- 0; probestofinterestDF$UPCgroup
# give them useful names
probestofinterestDF$UPCgroup [probestofinterestDF$UPCgroup ==1] <- "active"
probestofinterestDF$UPCgroup [probestofinterestDF$UPCgroup ==0] <- "inactive"

#  assign rownames (probeset) as unique column names and make it a factor
# (*critical* for ggplot2 scatter plot)
probestofinterestDF$sampleID <- row.names(probestofinterestDF)
probestofinterestDF$sampleID <- factor(probestofinterestDF$sampleID, levels = probestofinterestDF$sampleID)

# draw the plot
library(ggplot2)

# log values
lineintercept <- log(UPCthreshold)
probestofinterestDF.log <- probestofinterestDF
probestofinterestDF.log$probeset <- log(probestofinterestDF$probeset)
maxY <- max(probestofinterestDF.log$probeset, na.rm = TRUE) # close to log(UPC of 1)
minY <- min(probestofinterestDF.log$probeset, na.rm = TRUE) # close to log(UPC of 0)
lineintercept <- log(UPCthreshold)

rect <- data.frame(xmin=-Inf, xmax=Inf, ymin=(minY-1), ymax=log(0.1)) # draw a shaded area

UPCplot <- ggplot(probestofinterestDF.log, aes(x=sampleID, y=probeset)) + geom_point(size = 1, aes(color=factor(UPCgroup))) + geom_hline(aes(yintercept=lineintercept), colour="black", linetype="dashed") + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey20", alpha=0.2, inherit.aes = FALSE,linetype="blank") + theme_bw() + scale_color_manual(values=c("#51632c", "#c02026")) + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
# if you want to see sample text labels: +geom_text(aes(label=sampleID),hjust=0, vjust=0,size=0.2)

ggsave(filename="./scatter.plot.log.UPC.pdf", plot=UPCplot, width=15, height=5, dpi=1200)
# ~~~~~~~~~~~~~~~~~~~~~~~~ THIS IS THE END ~~~~~~~~~~~~~~~~~~~~~~~~ #
