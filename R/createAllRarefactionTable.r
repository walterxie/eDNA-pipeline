library(vegan)
library(Hmisc)


if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("isPlot")) isPlot = FALSE # by subplot
if(!exists("otuThr")) otuThr = 97
if(!exists("levels")) levels = rep(c("gamma","alpha","beta"),3)
if(!exists("qs")) qs = rep(0:2,each=3)

source("Modules/init.r")

######## get min size all sites for rdiversityTable each gene #######
getMinSizeAllSites <- function(communityMatrix) {
	minSample <- min(99999, min(rowSums(communityMatrix)))  		

	sampleSize <- minSample#floor(minSample/100) * 100 # make sure subsampling run  
	cat("Community matrix: min sample size per site =", minSample, ", take sample size per site =", sampleSize, "\n") 
    
    return (sampleSize)
}

cat("Intermediate data: create rarefaction diversity table using rmSingleton =", rmSingleton, ", isPlot =", isPlot, ", otuThr =", otuThr, "\n") 

n <- length(matrixNames)
# main
for (expId in 1:n) {	
  isP <- isPlot
  min2 <- rmSingleton
	if (expId == n) {
	  isP <- TRUE
	  min2 <- FALSE
	} 
  cat("Create rarefaction diversity table for", matrixNames[expId], 
      ", rmSingleton =", min2, ", isPlot =", isP, ", otuThr =", otuThr, ".\n")
  communityMatrix <- getCommunityMatrixT(expId, isP, min2)
	sampleSize <- getMinSizeAllSites(communityMatrix) 

	source("Modules/RarefactionTable.R", local=TRUE)

	# create rarefaction table file which is calculated by RarefactionTable.R	
	outputRFTable <- paste(workingPath, "data/", postfix(matrixNames[expId], isPlot, rmSingleton, sep="-"), "-rarefaction-table.csv", sep = "")
	writeRdiversityTable(communityMatrix, outputRFTable)
} 



