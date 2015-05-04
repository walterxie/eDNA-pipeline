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
getMinSizeAllSites <- function(expId, otuThr) {
    minSample = 99999

	communityMatrix <- init(expId, otuThr, stringBySubOrPlot)
	minSample <- min(minSample, min(rowSums(communityMatrix)))  		

	sampleSize <- minSample#floor(minSample/100) * 100 # make sure subsampling run  
	print(paste("otuThr =", otuThr, ", min sample size per site =", minSample, ", take sample size per site =", sampleSize)) 
    
    return (sampleSize)
}

# main
for (expId in 1:length(experiments)) {	
    print(paste("experiment = ", experiments[expId], ", matrixName = ", matrixNames[expId], sep=""))
    
	rdiversityOTUThreTableMean <- NULL
	
	# add "by-plot" to merge 2 nearby columns 
	matrixName <- paste(matrixNames[expId], stringBySubOrPlot, sep = "")
	
	# use same sample size for all OTU threshold, but diff for each gene
	minSample = 999999
	for(otuThr in 90:100) {
		minSample <- min(minSample, getMinSizeAllSites(expId, otuThr)) 
		sampleSize <- minSample
	}
	#sampleSize = 1100
	print(paste(experiments[expId], " takes sampleSize = ", sampleSize, sep=""))  


	for(otuThr in 90:100) {
		print(paste("Rarefraction at otuThr =", otuThr))
		communityMatrix <- init(expId, otuThr, stringBySubOrPlot)

		source("TurnoverDist.R", local=TRUE)
		# need to put on the top
		source("SampleCounts.R", local=TRUE)

		source("Diversities.R", local=TRUE)

        source("RarefractionDiversitiesTable.R", local=TRUE)

        if (otuThr == 97) {
			# create Rarefraction table file which is calculated by RarefractionDiversitiesTable.R
			outputRFTable <- paste(workingPath, experiments[expId], "/", matrixName, "-", otuThr, "-rarefraction-table.csv", sep = "")
			source("Rarefraction.R", local=TRUE) 	
		} 
	    
	    rfVsOTUReps <- 5
		if (is.null(rdiversityOTUThreTableMean)) {
			rdiversityOTUThreTableMean <- data.frame(row.names=paste(levels, qs, sep=""), check.names=FALSE)
			rdiversityOTUThreTableMean <- rdiversityTable(communityMatrix, sampleSize, rfVsOTUReps)[,1]
		} else {
			tmpRDT <- rdiversityTable(communityMatrix, sampleSize, rfVsOTUReps)[,1] 
			rdiversityOTUThreTableMean <- cbind(rdiversityOTUThreTableMean, tmpRDT)
		}			
	}

	colnames(rdiversityOTUThreTableMean) <- paste("otus.", 90:100, sep="")		

	outputRDT <- paste(workingPath, experiments[expId], "/", matrixName, "-otus-thre-table.csv", sep = "")
	write.csv(rdiversityOTUThreTableMean, outputRDT, row.names=FALSE, quote=FALSE)		
} 



