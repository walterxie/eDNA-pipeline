library(vegan)
library(vegetarian)
library(grid)
library(gridExtra)
library(data.table)
library(xtable)
library(tools)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(colorspace)
#library(Matrix)
library(reshape2)
library(ggdendro)
library(plyr)
library(cluster)
library(ape)

if(!exists("sourcePath")) stop("source path to initiate modules is missing !")
if(!exists("workingPath")) stop("working path containing data is missing !")

if(!exists("verbose")) verbose <- TRUE

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
# may contain empty rows cols
getCommunityMatrix <- function(expId, isPlot, min2) {
  n <- length(matrixNames)
  matrixName <- matrixNames[expId]
  
  # hard code to get Vegetation
  if (expId==n) {
    if (!isPlot)
      stop("Vegetation only has plot based community matrix !")
    inputCM <- file.path(workingPath, "data", "LBI_Trees_Saplings_SBA.csv")
  } else if (isPlot) {
    inputCM <- file.path(workingPath, "data", paste(matrixName, "by_plot.txt", sep="_"))
  } else {
    # e.g. data/16S.txt
    inputCM <- file.path(workingPath, "data", paste(matrixName, ".txt", sep=""))
  }
  
  communityMatrix <- readFile(inputCM)
  if(verbose) 
    cat("\nUpload community matrix : ", ncol(communityMatrix), "columns,", nrow(communityMatrix), "rows, from", inputCM, "\n") 
  
  if (min2) {
    singletons <- which(rowSums(communityMatrix)==1)
    if(verbose) 
      cat("Remove", length(singletons) ,"singletons from ", matrixName, " !\n")
    communityMatrix <- communityMatrix[-singletons,]
    rm(singletons)		
  }
  
  return(communityMatrix)
}

# remove empty rows cols
prepCommunityMatrix <- function(communityMatrix) {
  # filter column first to avoid empty rows after columns remvoed if vectorThr>0
  if(any(colSums(communityMatrix)== 0)) {
    communityMatrix <- rmVectorFromCM(communityMatrix, vectorThr=0, MARGIN=2)
    #stop("Invalid input: community matrix has empty samples !")
  }
  if(any(rowSums(communityMatrix)== 0)) {
    communityMatrix <- rmVectorFromCM(communityMatrix, vectorThr=0, MARGIN=1)
    #stop("Invalid input: community matrix has empty OTUs !")
  }
  
  return(communityMatrix)
}

# transposed CM for vegan, and remove empty rows cols 
# return(NULL) if nrow(taxaPaths) < minRow, default minRow=0
getCommunityMatrixT <- function(expId, isPlot, min2, taxa.group="all", minRow=0) {
  communityMatrix <- getCommunityMatrix(expId, isPlot, min2)
  
  n <- length(matrixNames)
  if (expId < n && taxa.group != "all") {
    ##### load data #####
    taxaPaths <- getTaxaPaths(expId, taxa.group)
    
    if (nrow(taxaPaths) < minRow) {
      cat("Warning: return NULL, because", nrow(taxaPaths), "row(s) classified as", taxa.group, "<", minRow, "threshold.\n")
      return(NULL)
    } else {
      # merge needs at least 2 cols 
      taxaAssgReads <- merge(communityMatrix, taxaPaths, by = "row.names")
      
      if (nrow(taxaAssgReads) < minRow) {
        cat("Warning: return NULL, because", nrow(taxaAssgReads), "row(s) in cm match", taxa.group, "<", minRow, "threshold.\n")
        return(NULL)
      }
      
      # move 1st col Row.names to row.names
      rownames(taxaAssgReads) <- taxaAssgReads[,"Row.names"]
      taxaAssgReads <- taxaAssgReads[,-1]
      # get CM
      taxaAssgReads <- taxaAssgReads[,1:ncol(communityMatrix)]
      
      cat("Merging", nrow(taxaAssgReads), "matched OTUs from", nrow(communityMatrix), "OTUs in matrix to", 
          nrow(taxaPaths), "taxa classification, taxa.group =", taxa.group, ", final ncol =", ncol(taxaAssgReads), ".\n")
      
      communityMatrix <- data.matrix(taxaAssgReads)
    }
  }
  
  communityMatrix <- prepCommunityMatrix(communityMatrix)
  
  communityMatrixT <- transposeCM(communityMatrix)
  
  return(communityMatrixT)
}
# rownames(communityMatrix) <- gsub("-(.*)|", "\\1", rownames(communityMatrix))

###### taxa assignment by reads #####
# "ARCHAEA", "BACTERIA", "CHROMISTA", "PROTOZOA", "FUNGI", "PLANTAE", "ANIMALIA", "EUKARYOTA", "PROKARYOTA", "PROTISTS"
# PROKARYOTA: all prokaryotes (Bacteria + Archaea)
# EUKARYOTA: all eukaryotes
# PROTISTS: "CHROMISTA|PROTOZOA" all micro-eukaryotes
getTaxaPaths <- function(expId, taxa.group="all", rank="kingdom") {
  inputTaxa <- file.path(workingPath, "taxonomy_tables", paste(matrixNames[expId], "taxonomy_table.txt", sep="_"))
  taxaPaths <- readTaxaFile(inputTaxa)	
  taxaPaths <- taxaPaths[order(rownames(taxaPaths)),]
  # make lower case to match ranks
  colnames(taxaPaths) <- tolower(colnames(taxaPaths))
  
  nTaxa=nrow(taxaPaths)
  ##### keep OTU rows contain given taxa belongTo ##### 
  if (taxa.group != "all") {
    # Exclude unassigned etc
    taxaPaths <- subset(taxaPaths, !(grepl("root|cellular organisms|No hits|Not assigned", taxaPaths[,"kingdom"])))  
    # (Only retain prokaryotes for 16S, eukaryotes for the other amplicons)  
    if (expId==1) {
      taxaPaths <- subset(taxaPaths, (grepl("BACTERIA|ARCHAEA", taxaPaths[,"kingdom"])))  
    } else {
      taxaPaths <- subset(taxaPaths, !(grepl("BACTERIA|ARCHAEA", taxaPaths[,"kingdom"])))
    }
    # Exclude probably bogus taxa
    taxaPaths <- subset(taxaPaths, !(grepl("Cnidaria|Brachiopoda|Echinodermata|Porifera", taxaPaths[,"phylum"])|
                                       grepl("Bivalvia|Teleostei|Elasmobranchii|Polyplacophora", taxaPaths[,"class"])|
grepl("Nudibranchia|Crocodylia|Serpentes|Testudines|Carnivora|Gymnophiona|Lagomorpha|Rodentia|Serpentes|Scorpiones", taxaPaths[,"order"])))
    
    if (taxa.group != "assigned") {
      if (toupper(taxa.group) == "PROKARYOTA" || toupper(taxa.group) == "EUKARYOTA") {
        taxaPaths <- subset(taxaPaths, grepl(taxa.group, taxaPaths[,"superkingdom"])) 
      } else if (toupper(taxa.group) == "PROTISTS") {
        taxaPaths <- subset(taxaPaths, grepl("CHROMISTA|PROTOZOA", taxaPaths[,"kingdom"])) 
      } else {
        taxaPaths <- subset(taxaPaths, grepl(taxa.group, taxaPaths[,rank])) 
      }
    }
  }
  
  if (nrow(taxaPaths) < 1)
    cat("Warning: cannot find", taxa.group, "from taxa path file", inputTaxa, "!")
  
  if(verbose && nrow(taxaPaths) < nTaxa) 
    cat("\nSelect", nrow(taxaPaths), "taxa classification, taxa.group =", taxa.group, ".\n") 
  
  return(taxaPaths)
}


# rankLevel: the taxa level in each bar
# groupLevel: used to assign colour for each group, and must higher than rankLevel
# taxa.group: keep OTU rows contain given taxa group, if "all", keep all
# return taxaAssgReads = CM + rankLevel + groupLevel
getTaxaAssgReads <- function(expId, isPlot, min2, rankLevel, groupLevel, taxa.group="all") {
  cat("Create taxonomy assignment for", matrixNames[expId], ".\n")
  
  ##### load data #####
  communityMatrix <- getCommunityMatrix(expId, isPlot, min2)
  communityMatrix <- prepCommunityMatrix(communityMatrix)
  
  communityMatrix <- communityMatrix[order(rownames(communityMatrix)),]
  communityMatrix <- communityMatrix[,order(colnames(communityMatrix))]
  
  taxaPaths <- getTaxaPaths(expId, taxa.group)
  
  ###### taxa assignment by reads #####
  if ( ! tolower(rankLevel) %in% tolower(colnames(taxaPaths)) ) 
    stop( paste("Column name", rankLevel, "not exist in taxa path file for ", matrixNames[expId]) )
  if (! tolower(groupLevel) %in% tolower(colnames(taxaPaths)) ) 
    stop( paste("Column name", groupLevel, "not exist in taxa path file for ", matrixNames[expId]) )
  
  colRankLevel <- which(tolower(colnames(taxaPaths))==tolower(rankLevel))
  colGroupLevel <- which(tolower(colnames(taxaPaths))==tolower(groupLevel))
  
  taxaAssgReads <- merge(communityMatrix, taxaPaths[,c(colRankLevel, colGroupLevel)], by = "row.names")
  
  taxaAssgReads[,rankLevel] <- gsub("(\\s\\[=.*\\])", "", taxaAssgReads[,rankLevel])
  taxaAssgReads[,groupLevel] <- gsub("(\\s\\[=.*\\])", "", taxaAssgReads[,groupLevel])
  
  cat("Merging:", nrow(taxaAssgReads), "OTUs are matched from", nrow(communityMatrix), "OTUs in matrix to", 
      nrow(taxaPaths), "taxa classification, taxa.group =", taxa.group, ".\n")
  
  return(taxaAssgReads)
}

getTaxaRef <- function() {
  tax_ref <- read.table(file.path(workingPath, "data", "New_taxonomy_from_PLOSONE_2015_fixed.txt"), 
                        header = TRUE, sep = "\t", quote = "", comment.char = "")
  # make lower case to match ranks
  colnames(tax_ref) <- tolower(colnames(tax_ref))
  
  # Remove quirks/questions in taxa ([= ...])
  tax_ref <- apply(tax_ref, 2, function(col) gsub("(\\s\\[=.*\\])", "", col))
}

###### Trees #####
getPhyloTree <- function(fNameStem) {
  inputT <- file.path(workingPath, "trees", paste(fNameStem, "tre", sep = "."))
  if (file.exists(inputT)) {
    cat("Load tree from", inputT, "\n") 
    tree <- read.tree(inputT)
    if(verbose) print(tree)
  } else {
    tree <- NULL
    cat("Warning: cannot find tree file:", inputT, "\n") 
  }
  tree
}

###### table to plot Phylo Rarefaction ##### 
getPhyloRareTable <- function(expId, isPlot, min2, taxa.group="assigned") {
  n <- length(matrixNames) 
  # hard code for Vegetation that only has plot and always keep singletons
  if (expId == n) {
    mid.name <- postfix("all", TRUE, FALSE, sep="-")
  } else {
    mid.name <- postfix(taxa.group, isPlot, min2, sep="-") 
  }
  
  inputT <- file.path(workingPath, "data", "pdrf", paste(matrixNames[expId], mid.name, "phylorare", "table.csv", sep="-"))
  if (file.exists(inputT)) {
    phylo.rare.df <- read.csv(file=inputT, head=TRUE, sep=",", row.names=1, check.names=FALSE)
    if(verbose) 
      cat("\nUpload phylo rarefaction table from", inputT, "\n") 
  } else {
    phylo.rare.df <- NULL
    cat("Warning: cannot find phylo rarefaction table", inputT, "\n") 
  }
  phylo.rare.df
}

###### table to plot Rarefaction ##### 
getRarefactionTableTaxa <- function(expId, isPlot, min2, taxa.group, div="alpha1") {
  pathFileStem <- file.path(workingPath, "data", "rf", paste(matrixNames[expId], 
                    postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-"))
  inputT <- paste(pathFileStem, "rare", div, "table.csv", sep = "-")
  if (file.exists(inputT)) {
    rare.df <- read.csv(file=inputT, head=TRUE, sep=",", row.names=1, check.names=FALSE)
    if(verbose) 
      cat("\nUpload rarefaction table per sample from", inputT, "\n") 
  } else {
    rare.df <- NULL
    cat("Warning: cannot find rarefaction table per sample", inputT, "\n") 
  }
  rare.df
}

getRarefactionTable <- function(expId, isPlot, min2) {
  n <- length(matrixNames) 
  matrixName <- matrixNames[expId]
  # hard code for Vegetation that only has plot and always keep singletons
  if (expId == n) {
    matrixName <- postfix(matrixName, TRUE, FALSE, sep="-")
  } else {
    matrixName <- postfix(matrixName, isPlot, min2, sep="-") 
  }
  
  inputRDT <- file.path(workingPath, "data", paste(matrixName, "rarefaction-table.csv", sep="-"))
  if(verbose) 
    cat("\nUpload rarefaction table : from", inputRDT, "\n") 
  
  rarefactionTable <- read.csv(file=inputRDT, head=TRUE, sep=",", row.names=paste(levels, qs, sep=""), check.names=FALSE)
}

###### dissimilarity matrix #####
# Dissimilarity matrix of paired samples
# diss.fun = "beta1-1", "jaccard", "horn.morisita"
getDissimilarityMatrix <- function(expId, isPlot, min2, diss.fun="beta1-1", taxa.group="all") {
  n <- length(matrixNames) 
  # hard code for Vegetation that only has plot and always keep singletons
  if (expId == n) {
    fname <- paste(matrixNames[expId], postfix("all", TRUE, FALSE, sep="-"), diss.fun, sep = "-")
  } else {
    fname <- paste(matrixNames[expId], postfix(taxa.group, isPlot, min2, sep="-"), diss.fun, sep = "-") 
  }
  
  inputB <- file.path(workingPath, "data", "dist", paste(fname, "csv", sep = "."))
  if(verbose) 
    cat("\nUpload", diss.fun, "matrix of", taxa.group, "taxa group(s) from", inputB, "\n") 
  
  diss.matrix <- readFile(file=inputB, sep=",")
  
  return(diss.matrix)
}

###### table to max remained diversity ##### 
getMaxRemainedDiversity <- function(lev.q, taxa.group="assigned") {
  inputT <- file.path(workingPath, "data", "maxrd", paste("max-div", lev.q, taxa.group,"table.csv", sep = "-"))
  if (file.exists(inputT)) {
    max.rd <- read.csv(file=inputT, head=TRUE, sep=",", row.names=1, check.names=FALSE)
    if(verbose) 
      cat("\nUpload max remained diversity table from", inputT, "\n") 
  } else {
    max.rd <- NULL
    cat("Warning: cannot find max remained diversity table", inputT, "\n") 
  }
  max.rd
}

getMaxRemainedDiversityRank <- function(lev.q, taxa.group="assigned") {
  inputT <- file.path(workingPath, "data", "maxrd", paste("max-div-rank", lev.q, taxa.group,"table.csv", sep = "-"))
  if (file.exists(inputT)) {
    max.rd <- read.csv(file=inputT, head=TRUE, sep=",", row.names=1, check.names=FALSE)
    if(verbose) 
      cat("\nUpload max remained diversity rank table from", inputT, "\n") 
  } else {
    max.rd <- NULL
    cat("Warning: cannot find max remained diversity rank table", inputT, "\n") 
  }
  max.rd
}

######## meta data of samples #######
getSampleMetaData <- function(isPlot) {
  if (isPlot) {
    inputCM <- file.path(workingPath, "Env_data", "LBI_all_env_data_by_plot.txt")
  } else {
    # e.g. data/16S.txt
    inputCM <- file.path(workingPath, "Env_data", "LBI_all_env_data_by_subplot.txt")
  }
  if(verbose) 
    cat("\nUpload enviornmental data from", inputCM, "\n") 
  env <- readFile(inputCM)
}


######## elevations #######
getElevPlotDist <- function(plot.names, env.byplot) { 
  colElev = 1
  # case insensitive
  matched.id <- match(tolower(plot.names), tolower(rownames(env.byplot)))
  matched.id <- matched.id[!is.na(matched.id)]
  # match 
  env.plot.match <- env.plot[matched.id, ]
  
  cat("Find", nrow(env.plot.match), "plots having elevations, community matrix has", 
      length(plot.names), "plots, meta-data file has", nrow(env.plot), "plots.\n")
  
  return(dist(env.plot.match[,colElev]))
}


