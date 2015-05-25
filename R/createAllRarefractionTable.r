library(vegan)
library(Hmisc)

replicate = 1
threshold = 150
levels = rep(c("gamma","alpha","beta"),3)
qs = rep(0:2,each=3)

# change config below
sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
setwd(sourcePath)

workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
experiments <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name   
matrixNames <-  experiments # 

stringBySubOrPlot <- "-by-subplot"

source("init.R", local=TRUE)

# get min size all sites for rdiversityTable each gene
getMinSizeAllSites <- function(expId, otuThre) {
    minSample = 99999

	communityMatrix <- init(expId, otuThre, stringBySubOrPlot)
	minSample <- min(minSample, min(rowSums(communityMatrix)))  		

	sampleSize <- minSample#floor(minSample/100) * 100 # make sure subsampling run  
	print(paste("otuThre =", otuThre, ", min sample size per site =", minSample, ", take sample size per site =", sampleSize)) 
    
    return (sampleSize)
}

otuThre=97

# main
for (expId in 1:length(experiments)) {	
    print(paste("experiment = ", experiments[expId], ", matrixName = ", matrixNames[expId], sep=""))
    	
	# add "by-plot" to merge 2 nearby columns 
	matrixName <- paste(matrixNames[expId], stringBySubOrPlot, sep = "")
	
	# use same sample size for all OTU threshold, but diff for each gene
	sampleSize <- getMinSizeAllSites(expId, otuThre) 
	
	#sampleSize = 1100
	print(paste(experiments[expId], " takes sampleSize = ", sampleSize, sep=""))  

	print(paste("Rarefraction at otuThre =", otuThre))
	communityMatrix <- init(expId, otuThre, stringBySubOrPlot)

	source("TurnoverDist.R", local=TRUE)
	# need to put on the top
	source("SampleCounts.R", local=TRUE)

	source("Diversities.R", local=TRUE)

	source("RarefractionDiversitiesTable.R", local=TRUE)

	# create Rarefraction table file which is calculated by RarefractionDiversitiesTable.R
	outputRFTable <- paste(workingPath, experiments[expId], "/", matrixName, "-", otuThre, "-rarefraction-table.csv", sep = "")
	source("Rarefraction.R", local=TRUE) 
} 



