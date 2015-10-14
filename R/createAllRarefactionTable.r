library(vegan)
library(Hmisc)

# change config below
#sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
#setwd(sourcePath)
#workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
#matrixNames <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name   
#levels = rep(c("gamma","alpha","beta"),3)
#qs = rep(0:2,each=3)

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

cat("Analysis: create rarefaction diversity table using rmSingleton =", rmSingleton, ", isPlot =", isPlot, ", otuThr =", otuThr, "\n") 

# main
for (expId in 1:length(matrixNames)) {	
  cat("Calculate rarefaction table for", matrixNames[expId], ".\n")
  
	communityMatrix <- getCommunityMatrixT(expId, isPlot, rmSingleton)
	sampleSize <- getMinSizeAllSites(communityMatrix) 

	source("Modules/RarefactionTable.R", local=TRUE)

	# create rarefaction table file which is calculated by RarefactionTable.R	
	outputRFTable <- paste(workingPath, "data/", postfix(matrixNames[expId], isPlot, rmSingleton, sep="-"), "-rarefaction-table.csv", sep = "")
	writeRdiversityTable(communityMatrix, outputRFTable)
} 



