# create SBA matrix from DBH  
# http://oregonstate.edu/instruct/bot440/wilsomar/Content/HTM-trees.htm
#
# Author: Walter Xie
# Accessed on 29 Oct 2015

library(reshape)

# change path here
workingPath <- "~/WorkSpace/eDNA-pipeline/Vegetation"
setwd(workingPath)

###### load DBH, file format is: 
#Plot No	Plot Name	Species	Stem Count	DBH
#1	Plot01	BRAREP	1	1.3
#1	Plot01	BRAREP	1	4
dbhFile <- file.path(workingPath, "LBI_Tree_DBH_Values_All_Trees_Saplings.csv")
dbh <- read.csv(dbhFile, head=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
cat("\nUpload DBH file: ", nrow(dbh), "rows", ncol(dbh), "columns from", dbhFile, "\n")

col_names <- c("PlotName","Species","StemCount","DBH")
colnames(dbh)

if ( ! all(col_names %in% colnames(dbh)) ) 
  stop("Invalid file format !")

###### calculate SBA
A <- 400 # m2
dbh$BA <- pi*(dbh$DBH/2)^2 # cm2
dbh$BA.sum <- dbh$BA * dbh$StemCount

sba <- cast(dbh, Species~PlotName, sum)

# move 1st col to row names
rownames(sba) <- sba[,"Species"]
sba <- sba[,-which(colnames(sba) %in% "Species")]

sba <- sba / A # cm2/m2 or m2/ha

###### write matrix to file
sbaFile <- file.path(workingPath, "LBI_Trees_Saplings_SBA.csv")
cat("\nWrite SBA matrix having", nrow(sba), "species", ncol(sba), "plots to file", sbaFile, ".\n")
write.csv(sba, sbaFile, quote=FALSE)

