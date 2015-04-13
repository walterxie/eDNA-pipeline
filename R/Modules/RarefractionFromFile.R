# returns a table of diversities for the given community matrix
# and sample size to rarefy to
# requires the d() function from vegetarian library
# requires source("Modules/TurnoverDist.R"), 
# source("Modules/SampleCounts.R"), source("Modules/Diversities.R")

levels = rep(c("gamma","alpha","beta"),3)
qs = rep(0:2,each=3)

if(!exists("inputRDT")) stop("input path for rarefraction table is missing !")
# 9 rows of Jost diversities, columns are mean, sterrD   
rdiversityTableMean <- read.csv(file=inputRDT, head=TRUE, sep=",", row.names=paste(levels, qs, sep=""), check.names=FALSE)

print(paste("input table of diversities vs. sample size per site from file : ", inputRDT, sep="")) 

# remove prefix size.
sampleSizesSeq <- gsub("^.*?\\.","",colnames(rdiversityTableMean)) 
	    
