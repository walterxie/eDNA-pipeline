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

otuThr = 97

stringBySubOrPlot <- "-by-subplot"

source("Modules/init.R", local=TRUE)

# main
for (expId in 1:length(matrixNames)) {	
    print(paste("experiment = ", matrixNames[expId], ", matrixName = ", matrixNames[expId], sep=""))
    	
	# add "by-plot" to merge 2 nearby columns 
	matrixName <- paste(matrixNames[expId], stringBySubOrPlot, sep = "")
	if (rmSingleton) 
		matrixName <- paste("nonsingleton", matrixName, sep = "-")
			
	# use same sample size for all OTU threshold, but diff for each gene
	sampleSize <- getMinSizeAllSites(expId, otuThre) 
	
	#sampleSize = 1100
	print(paste(matrixNames[expId], " takes sampleSize = ", sampleSize, sep=""))  

	print(paste("Rarefaction at otuThre =", otuThre))
	communityMatrix <- init(expId, otuThre, stringBySubOrPlot)

	source("Modules/TurnoverDist.R", local=TRUE)
	# need to put on the top
	source("Modules/SampleCounts.R", local=TRUE)

	source("Modules/RarefactionDiversitiesTable.R", local=TRUE)

	# create rarefaction table file which is calculated by RarefactionDiversitiesTable.R	
	outputRFTable <- paste(workingPath, matrixNames[expId], "/", matrixName, "-", otuThre, "-rarefaction-table.csv", sep = "")
	source("Modules/Rarefaction.R", local=TRUE) 
} 



