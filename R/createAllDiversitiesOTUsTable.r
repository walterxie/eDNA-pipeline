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

stringBySubOrPlot <- "-by-subplot"

source("Modules/init.R", local=TRUE)

# main
for (expId in 1:length(matrixNames)) {	
    print(paste("experiment = ", matrixNames[expId], ", matrixName = ", matrixNames[expId], sep=""))
    
	rdiversityOTUThreTableMean <- NULL
	
	# add "by-plot" to merge 2 nearby columns 
	matrixName <- paste(matrixNames[expId], stringBySubOrPlot, sep = "")
	if (rmSingleton) 
		matrixName <- paste("nonsingleton", matrixName, sep = "-")
	
	# use same sample size for all OTU threshold, but diff for each gene
	minSample = 999999
	for(otuThr in 90:100) {
		minSample <- min(minSample, getMinSizeAllSites(expId, otuThr)) 
		sampleSize <- minSample
	}
	#sampleSize = 1100
	print(paste(matrixNames[expId], " takes sampleSize = ", sampleSize, sep=""))  


	for(otuThr in 90:100) {
		print(paste("Rarefaction at otuThr =", otuThr))
		communityMatrix <- init(expId, otuThr, stringBySubOrPlot)

		source("Modules/TurnoverDist.R", local=TRUE)
		# need to put on the top
		source("Modules/SampleCounts.R", local=TRUE)

        source("Modules/RarefactionDiversitiesTable.R", local=TRUE)

        if (otuThr == 97) {
			# create rarefaction table file which is calculated by RarefactionDiversitiesTable.R
			outputRFTable <- paste(workingPath, matrixNames[expId], "/", matrixName, "-", otuThr, "-rarefaction-table.csv", sep = "")
			source("Modules/Rarefaction.R", local=TRUE) 	
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

	outputRDT <- paste(workingPath, matrixNames[expId], "/", matrixName, "-otus-thre-table.csv", sep = "")
	write.csv(rdiversityOTUThreTableMean, outputRDT, row.names=FALSE, quote=FALSE)		
} 



