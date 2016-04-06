
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("taxa.group")) taxa.group="assigned"
if(!exists("isPlot")) isPlot = TRUE # by subplot


n <- length(matrixNames) 

source("Modules/init.r")
source("Modules/MaximizeDiversity.R")

filePath <- file.path(workingPath, "data", "maxrd")
mkdir(filePath) 

# no singlton and plot only
cat("Intermediate data: create max remained diversity table",
    ", rmSingleton =", TRUE, ", isPlot =", TRUE, ", otuThr =", otuThr, "\n") 

q=1
for (lev in c("gamma", "beta")) {
  if (lev == "beta") {
    q=1
  }
  
  # ranks
  allRanks <- NULL 
  # max diversity
  allMaxDiv <- NULL
  
  ######## eDNA #######
  for (expId in 1:(n-1)) {	
    communityMatrix <- getCommunityMatrixT(expId, isPlot=TRUE, min2=TRUE, taxa.group, minRow=600)
    
    if (is.null(communityMatrix)) {
      cat("\nWarning: OTUs < 600, skip", taxa.group, "subset from", matrixNames[expId], ".\n")
      next #return(FALSE)
    }
    
    # 3 columns: rank, diversity, site
    maxDiv <- getMaxRemainedDiversity(communityMatrix, level, q)
    
    if (is.null(allRanks)) {		
      allRanks <- data.frame(row.names=rownames(maxDiv), stringsAsFactors = FALSE)
      allMaxDiv <- data.frame(row.names=rownames(maxDiv), stringsAsFactors = FALSE)		
    } else {
      if(!all(rownames(allRanks) == rownames(maxDiv))) 
        stop(paste("Find different plots ", rownames(allRanks), " != ", rownames(maxDiv), sep=""))
    }
    
    allRanks[,matrixNames[expId]] <- as.numeric(maxDiv[,1]) 
    allMaxDiv[,matrixNames[expId]] <- as.numeric(maxDiv[,2])
  }  # END for expId
  
  if (!is.null(allRanks)) {
    ######## rank table #######
    allRanks$Mean <- rowMeans(allRanks, na.rm=T)
    allRanks$Sd <- apply(allRanks[,1:(ncol(allRanks)-1)], 1, sd)
    
    allRanks$Mean <- round(allRanks$Mean, 2)
    allRanks$Sd <- round(allRanks$Sd, 2)
    
    outputFile <- file.path(filePath, paste("max-div-rank", lev, q, taxa.group,"table.csv", sep = "-"))
    write.csv(allRanks, outputFile, quote=FALSE)
    outputFile <- file.path(filePath, paste("max-div", lev, q, taxa.group,"table.csv", sep = "-"))
    write.csv(allMaxDiv, outputFile, quote=FALSE)
    
    cat("Write max remained diversity table", nrow(allMaxDiv), "rows, ", ncol(allMaxDiv), "columns. \n") 
  }
}

