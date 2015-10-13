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
if(!exists("levels")) stop("levels of Jost diversity are missing !")
if(!exists("qs")) stop("qs of Jost diversity are missing !")

if(!exists("otuThr")) otuThr = 97

isPlot <- FALSE # by subplot

source("Modules/init.r")

######## get min size all sites for rdiversityTable each gene #######
getMinSizeAllSites <- function(communityMatrix) {
	minSample <- min(99999, min(rowSums(communityMatrix)))  		

	sampleSize <- minSample#floor(minSample/100) * 100 # make sure subsampling run  
	print(paste("Community matrix: min sample size per site =", minSample, ", take sample size per site =", sampleSize)) 
    
    return (sampleSize)
}

# main
for (expId in 1:length(matrixNames)) {	
    print(paste("experiment = ", matrixNames[expId], ", matrixName = ", matrixNames[expId], sep=""))
    	
	matrixName <- matrixNames[expId]
	if (rmSingleton) 
		matrixName <- paste(matrixName, "min2", sep = "-")
			
	communityMatrix <- getCommunityMatrixT(expId, isPlot)
	sampleSize <- getMinSizeAllSites(communityMatrix) 

	source("Modules/RarefactionTable.R", local=TRUE)

	# create rarefaction table file which is calculated by RarefactionTable.R	
	outputRFTable <- paste(workingPath, "data/", matrixName, "-rarefaction-table.csv", sep = "")
	writeRdiversityTable(communityMatrix, outputRFTable)
} 



