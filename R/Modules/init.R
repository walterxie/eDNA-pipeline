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
#library(Matrix)

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

# get plot names from subplots vector separated by sep
getPlot <- function(subplots, sep="-") {
  sapply(strsplit(as.character(subplots), sep), "[[", 1)
}

######## load community matrix #######
# expId = 1:6
# isPlot determines to use which matrix file, by subplot or plot 
# min2 = rmSingleton, whether remove all singletons
getCommunityMatrixT <- function(expId, isPlot, min2) {
  if(!exists("verbose")) verbose <- FALSE 
  
  n <- length(matrixNames)
  matrixName <- matrixNames[expId]
  
  # hard code to get Vegetation
  if (expId==n) {
    if (!isPlot)
      stop("Vegetation only has plot based community matrix !")
    inputCM <- paste(workingPath, "data/trees_saplings_by_plot.txt", sep="")
  } else if (isPlot) {
    inputCM <- paste(workingPath, "data/", matrixName, "_by_plot.txt", sep="")
  } else {
    # e.g. data/16S.txt
    inputCM <- paste(workingPath, "data/", matrixName, ".txt", sep="")
  }
  
  communityMatrix <- readFile(inputCM)
  if(verbose) 
    cat("\nUpload community matrix : ", ncol(communityMatrix), "columns,", nrow(communityMatrix), "rows, from", inputCM, "\n") 
  
  # filter column first to avoid empty rows after columns remvoed
  if(any(colSums(communityMatrix)== 0))
    communityMatrix <- rmVectorFromCM(communityMatrix, vectorThr=0, MARGIN=2)
  #stop("Invalid input: community matrix has empty samples !")
  if(any(rowSums(communityMatrix)== 0))
    communityMatrix <- rmVectorFromCM(communityMatrix, vectorThr=0, MARGIN=1)
  #stop("Invalid input: community matrix has empty OTUs !")
  
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
# rownames(communityMatrix) <- gsub("-(.*)|", "\\1", rownames(communityMatrix))

###### table to plot Rarefaction ##### 
getRarefactionTable <- function(expId, isPlot, min2) {
  n <- length(matrixNames) 
  matrixName <- matrixNames[expId]
  # hard code for Vegetation that only has plot and always keep singletons
  if (expId == n) {
    matrixName <- postfix(matrixName, TRUE, FALSE, sep="-")
  } else {
    matrixName <- postfix(matrixName, isPlot, min2, sep="-") 
  }
  
  inputRDT <- paste(workingPath, "data/", matrixName, "-rarefaction-table.csv", sep = "")   
  if(verbose) 
    cat("\nUpload rarefaction table : from", inputRDT, "\n") 
  
  rarefactionTable <- read.csv(file=inputRDT, head=TRUE, sep=",", row.names=paste(levels, qs, sep=""), check.names=FALSE)
}

###### taxa assignment by reads #####
unclassTaxa <- c("Not assigned", "No hits", "cellular organisms", "root")
# rankLevel: the taxa level in each bar
# groupLevel: used to assign colour for each group, and must higher than rankLevel
# belongTo: keep OTU rows contain given taxa belongTo, if NULL, keep all
# return CM + rankLevel + groupLevel
getTaxaAssgReads <- function(expId, min2, rankLevel, groupLevel, belongTo) {
  if(missing(belongTo)) belongTo<-NULL
  
  cat("Create taxonomy assignment for", matrixNames[expId], ".\n")
  
  ##### load data #####
  communityMatrix <- getCommunityMatrixT(expId, TRUE, min2)  
  # rotate to make rows to be OTUs to match taxa file
  communityMatrix <- transposeCM(communityMatrix) 
  
  communityMatrix <- communityMatrix[order(rownames(communityMatrix)),]
  communityMatrix <- communityMatrix[,order(colnames(communityMatrix))]
  
  inputTaxa <- paste(workingPath, "taxonomy_tables/", matrixNames[expId], "_taxonomy_table.txt", sep="")
  taxaPaths <- readTaxaFile(inputTaxa)	
  taxaPaths <- taxaPaths[order(rownames(taxaPaths)),]
  # make lower case to match ranks
  colnames(taxaPaths) <- tolower(colnames(taxaPaths))
  # define unclassified
  taxaPaths[taxaPaths==unclassTaxa[1] | taxaPaths==unclassTaxa[2] | taxaPaths==unclassTaxa[3] | is.na(taxaPaths)] <- "unclassified"
  
  if ( ! tolower(rankLevel) %in% tolower(colnames(taxaPaths)) ) 
    stop( paste("Column name", rankLevel, "not exist in taxa path file for ", matrixNames[expId]) )
  if (! tolower(groupLevel) %in% tolower(colnames(taxaPaths)) ) 
    stop( paste("Column name", groupLevel, "not exist in taxa path file for ", matrixNames[expId]) )
  
  ##### keep OTU rows contain given taxa belongTo ##### 
  if (!is.null(belongTo)) {
    taxaPaths <- taxaPaths[which(grepl(belongTo, taxaPaths[,1])),] # taxaPaths[,1] is taxa path separated by ;
    if (belongTo == "Metazoa")  # non Arthropoda
      taxaPaths <- taxaPaths[-which(grepl("Arthropoda", taxaPaths[,1])),]		
  }
  
  ###### taxa assignment by reads #####
  colRankLevel <- which(tolower(colnames(taxaPaths))==tolower(rankLevel))
  colGroupLevel <- which(tolower(colnames(taxaPaths))==tolower(groupLevel))
  
  taxaAssgReads <- merge(communityMatrix, taxaPaths[,c(colRankLevel, colGroupLevel)], by = "row.names")
  
  cat("Merging:", nrow(taxaAssgReads), "OTUs are matched from", nrow(communityMatrix), "in matrix to", 
      nrow(taxaPaths), "taxa classification.\n")
  
  return(taxaAssgReads)
}

######## beta1-1 #######
getBeta1Minus1 <- function(expId, isPlot, min2) {
  n <- length(matrixNames) 
  matrixName <- matrixNames[expId]
  # hard code for Vegetation that only has plot and always keep singletons
  if (expId == n) {
    matrixName <- postfix(matrixName, TRUE, FALSE, sep="-")
  } else {
    matrixName <- postfix(matrixName, isPlot, min2, sep="-") 
  }
  
  inputB <- paste(workingPath, "data/", matrixName, "-beta1-1.csv", sep = "")
  if(verbose) 
    cat("\nUpload beta1-1 matrix : from", inputB, "\n") 
  
  beta1_1 <- readFile(file=inputB, sep=",")
  
  return(beta1_1)
}


######## meta data of samples #######
getSampleMetaData <- function(isPlot) {
  if (isPlot) {
    inputCM <- paste(workingPath, "Env_data/LBI_all_env_data_by_plot.txt", sep="")
  } else {
    # e.g. data/16S.txt
    inputCM <- paste(workingPath, "Env_data/LBI_all_env_data_by_plot.txt", sep="")
  }
  
  env <- readFile(inputCM)
}


######## elevations #######
getElevPlotDist <- function(plot.names, env.plot) { 
  colElev = 1
  
  matched.id <- match(plot.names, rownames(env.plot))
  matched.id <- matched.id[!is.na(matched.id)]
  # match 
  env.plot.match <- env.plot[matched.id, ]

  cat("Find", nrow(env.plot.match), "plots having elevations, community matrix has", 
      length(plot.names), "plots, meta-data file has", nrow(env.plot), "plots.\n")

  return(dist(env.plot.match[,colElev]))
}

getElevSample <- function(cmRowNames) { 
  colPlot = 1
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
      elevTmp[i,colPlot] <- unlist(strsplit(as.character(elevTmp[i,colPlot]), "-", fixed = TRUE))[1] # hard code for names Plot1-B
      elevTmp[i+1,colPlot] <- unlist(strsplit(as.character(elevTmp[i+1,colPlot]), "-", fixed = TRUE))[1] 
      if(elevTmp[i,colPlot]  != elevTmp[i+1,colPlot] ) stop("Two subplots within the same plot should be next each other in the meta data file !") 
    }
    
    # map subplots meta data to plots
    elevTmp <- elevTmp[-seq(1,nrow(elevTmp), by=2),]		
  } 
  
  # missing data from traditional method
  if (sample_count == 8) {
    print(paste("Remove sample ", elevTmp[c(7,8),colPlot], " from meta data file due to miss data."))
    elevTmp <- elevTmp[c(-7,-8),]
  } 
  
  # remove not mapped samples from tmpEnvData
  if ( all(cmRowNames != elevTmp[,colPlot]) ) {
    stop(paste("Incorrect plots mapping between CM and meta data :  cm = ", cmRowNames, ", meta data =", elevTmp[,colPlot]))
  } 
  
  elevPlot <- data.frame(row.names=cmRowNames)
  elevPlot$Elevation.m <- elevTmp[,colElev]
  return(elevPlot) # 2nd col Elevation.m
}

