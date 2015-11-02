# create Stem Counts matrix from DBH  
# http://oregonstate.edu/instruct/bot440/wilsomar/Content/HTM-trees.htm
#
# Author: Walter Xie
# Accessed on 29 Oct 2015

library(reshape)

# change path here
workingPath <- "~/WorkSpace/eDNA-pipeline/Vegetation"
setwd(workingPath)

###### load DBH, file format is: 
#PlotNo	PlotName	Species	StemCount	DBH
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
dbh.threshold <- 10 # >= 10 cm
dbh <- dbh[dbh$DBH >= dbh.threshold,]
# make last col StemCount
dbh <- dbh[,-ncol(dbh)]

sc <- cast(dbh, Species~PlotName, sum)

# move 1st col to row names
rownames(sc) <- sc[,"Species"]
sc <- sc[,-which(colnames(sc) %in% "Species")]

###### write matrix to file
scFile <- file.path(workingPath, "LBI_Tree_Stem_Counts_10cm_DBH.csv")
cat("\nWrite stem counts matrix (DBH >= 10cm) having", nrow(sc), "species", ncol(sc), "plots to file", scFile, ".\n")
write.csv(sc, scFile, quote=FALSE)

