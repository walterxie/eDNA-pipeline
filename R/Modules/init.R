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

if(!exists("sourcePath")) stop("source path to initiate modules is missing !")
if(!exists("workingPath")) stop("working path containing data is missing !")

# utils, run them before analysis
source("Modules/IO.R")
source("Modules/Utils.R")
source("Modules/UtilsCM.R")
# fundmental functions
source("Modules/JostDiversities.R")

######## load community matrix #######
# expId = 1:6
init <- function(expId, otuThr, stringBySubOrPlot) {    
	matrixName <- matrixNames[expId]
	
	# e.g. otus97/16S_97_cm.csv
	inputCM <- paste(workingPath, matrixName, "/otus", otuThr, "/", matrixName, "-", otuThr, ".csv", sep="")
	
	#print(paste("expId = ", expId, ", matrixName = ", matrixName, ", inputCM = ", inputCM, sep=""))  
	
	matrixName <- paste(matrixName, stringBySubOrPlot, sep = "")

	source("Modules/InputBySubplot.R", local=TRUE)
	
	source("Modules/CommunityTransposeFilter.R", local=TRUE)
	
	if (rmSingleton) {
		singletons <- which(colSums(communityMatrix)==1)
		print(paste("Remove", length(singletons) ,"singletons from ", matrixName, " ! "))
		communityMatrix <- communityMatrix[,-singletons]
		rm(singletons)		
	}
	
	return (communityMatrix)
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

