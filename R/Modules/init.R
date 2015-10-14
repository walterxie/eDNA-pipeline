library(vegan)
library(vegetarian)
library(grid)
library(gridExtra)
library(data.table)
library(xtable)
library(tools)
library(tidyr)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(colorspace)
library(Matrix)

if(!exists("sourcePath")) stop("source path to initiate modules is missing !")
if(!exists("workingPath")) stop("working path containing data is missing !")

# utils, run them before analysis
source("Modules/IO.R")
source("Modules/Utils.R")
source("Modules/UtilsCM.R")
# fundmental functions
source("Modules/Diversities.R")

# most abundant 150 OTUs
threshold = 150

# add postfix for figure name or table label
postfix <- function(name, isPlot, min2, sep) {
	if (isPlot) 
	  name <- paste(name, "byplot", sep = sep) 
	if (min2) 
	  name <- paste(name, "min2", sep = sep)
	
	return(name)
}

######## load community matrix #######
# expId = 1:6
# isPlot determines to use which matrix file, by subplot or plot 
# min2 = rmSingleton, whether remove all singletons
getCommunityMatrixT <- function(expId, isPlot, min2) {
	if(!exists("verbose")) verbose <- FALSE 
    
	matrixName <- matrixNames[expId]
	
	if (isPlot) {
		inputCM <- paste(workingPath, "data/", matrixName, "_by_plot.txt", sep="")
	} else {
		# e.g. data/16S.txt
		inputCM <- paste(workingPath, "data/", matrixName, ".txt", sep="")
	}
	
	communityMatrix <- readFile(inputCM)
	if(verbose) 
		cat("\nUpload community matrix : ", ncol(communityMatrix), "columns,", nrow(communityMatrix), "rows, from", inputCM, "\n") 
	
	communityMatrixT <- transposeCM(communityMatrix)
	
	if (min2) {
		singletons <- which(colSums(communityMatrixT)==1)
		if(verbose) 
			cat("Remove", length(singletons) ,"singletons from ", matrixName, " !\n")
		communityMatrixT <- communityMatrixT[,-singletons]
		rm(singletons)		
	}
	
	return(communityMatrixT)
}

# table to plot Rarefaction
getRarefactionTable <- function(expId, isPlot, min2) {
    matrixName <- matrixNames[expId]
	if (isPlot) 
	  matrixName <- paste(matrixName, "byplot", sep = "-") 
	if (rmSingleton) 
	  matrixName <- paste(matrixName, "min2", sep = "-")  
			
	inputRDT <- paste(workingPath, "data/", matrixName, "-rarefaction-table.csv", sep = "")    
    rarefactionTable <- read.csv(file=inputRDT, head=TRUE, sep=",", row.names=paste(levels, qs, sep=""), check.names=FALSE)
}


######## elevations #######
# plot_elevations.txt 
#SampleName	Elevation.m		slope.degree	aspect.degree
#Plot1-B	50	18	270
#Plot1-L	50	18	270
#Plot2-C	90	30	60
#Plot2-M	90	30	60
#Plot3-H	160	36	120
#Plot3-N	160	36	120
getElevSampleDist <- function(cmRowNames) { 
    colSample = 1
    colElev = 2
    
    sample_count <- length(cmRowNames)
	print(paste("sample_count = ", sample_count, sep=""))

	inputElevation <- paste(workingPath, "data/plot_elevations.txt", sep="")	
	elev <- read.table(inputElevation, header=TRUE, stringsAsFactors=FALSE)

    # keep elev for easy debug
	elevTmp <- elev
	
   # hard code to fit in given data
    if (sample_count <= 10) {
		# 2 subplots within the plot should be next each other 
		for (i in seq(1,nrow(elevTmp),by=2)) {
			elevTmp[i,colSample] <- unlist(strsplit(as.character(elevTmp[i,colSample]), "-", fixed = TRUE))[1] # hard code for names Plot1-B
			elevTmp[i+1,colSample] <- unlist(strsplit(as.character(elevTmp[i+1,colSample]), "-", fixed = TRUE))[1] 
			if(elevTmp[i,colSample]  != elevTmp[i+1,colSample] ) stop("Two subplots within the same plot should be next each other in the meta data file !") 
		}
	  		
	  	# map subplots meta data to plots
		elevTmp <- elevTmp[-seq(1,nrow(elevTmp), by=2),]		
	} 
	
	# missing data from traditional method
	if (sample_count == 8) {
	  print(paste("Remove sample ", elevTmp[c(7,8),colSample], " from meta data file due to miss data."))
	  elevTmp <- elevTmp[c(-7,-8),]
	} 

    # remove not mapped samples from tmpEnvData
	if ( all(cmRowNames != elevTmp[,colSample]) ) {
		stop(paste("Incorrect plots mapping between CM and meta data :  cm = ", cmRowNames, ", meta data =", elevTmp[,colSample]))
	} 

	elevPlotDist <- dist(elevTmp[,colElev]) # 2nd col Elevation.m
}
######## elevations #######
getElevSample <- function(cmRowNames) { 
    colSample = 1
    colElev = 2
    
    sample_count <- length(cmRowNames)
	print(paste("sample_count = ", sample_count, sep=""))

	inputElevation <- paste(workingPath, "data/plot_elevations.txt", sep="")	
	elev <- read.table(inputElevation, header=TRUE, stringsAsFactors=FALSE)

    # keep elev for easy debug
	elevTmp <- elev
	
   # hard code to fit in given data
    if (sample_count <= 10) {
		# 2 subplots within the plot should be next each other 
		for (i in seq(1,nrow(elevTmp),by=2)) {
			elevTmp[i,colSample] <- unlist(strsplit(as.character(elevTmp[i,colSample]), "-", fixed = TRUE))[1] # hard code for names Plot1-B
			elevTmp[i+1,colSample] <- unlist(strsplit(as.character(elevTmp[i+1,colSample]), "-", fixed = TRUE))[1] 
			if(elevTmp[i,colSample]  != elevTmp[i+1,colSample] ) stop("Two subplots within the same plot should be next each other in the meta data file !") 
		}
	  		
	  	# map subplots meta data to plots
		elevTmp <- elevTmp[-seq(1,nrow(elevTmp), by=2),]		
	} 
	
	# missing data from traditional method
	if (sample_count == 8) {
	  print(paste("Remove sample ", elevTmp[c(7,8),colSample], " from meta data file due to miss data."))
	  elevTmp <- elevTmp[c(-7,-8),]
	} 

    # remove not mapped samples from tmpEnvData
	if ( all(cmRowNames != elevTmp[,colSample]) ) {
		stop(paste("Incorrect plots mapping between CM and meta data :  cm = ", cmRowNames, ", meta data =", elevTmp[,colSample]))
	} 

    elevPlot <- data.frame(row.names=cmRowNames)
    elevPlot$Elevation.m <- elevTmp[,colElev]
	return(elevPlot) # 2nd col Elevation.m
}

