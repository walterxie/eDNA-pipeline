# "-by-plot" trigger merge 2 subplots columns

# change config below
sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
setwd(sourcePath)

workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
experiments <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name
matrixNames <-  experiments

matrixNamesNo454 <-  c("seedlings","trees","invertebrates","birds") # need expId correct for "birds","seedlings" 

n <- length(matrixNames) 
m <- length(matrixNamesNo454)

levels = rep(c("gamma","alpha","beta"),3)
qs = rep(0:2,each=3)

otuThr = 97

source("init.R", local=TRUE)

######## 454 #######	

allStatTable <- NULL
allStatTable$stat <- c("unique reads","OTUs","alpha","effective alpha")

for (expId in 1:n) {	
    communityMatrix <- init(expId, otuThr, "-by-plot")
    
	source("TurnoverDist.R", local=TRUE)
	# need to put on the top
	source("SampleCounts.R", local=TRUE)

	source("Diversities.R", local=TRUE)

	print(paste(matrixNames[expId], " (", otuThr, "%) : OTUs = ", ncol(communityMatrix), ", total reads = ", sum(communityMatrix), "; min sample per plot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))
	print(paste("alpha = ", diversity_table[2,1], "; effective alpha = ", diversity_table[2,2], sep=""))
	print(diversity_table)
	
	allStatTable$tmp <- c(sum(communityMatrix),ncol(communityMatrix),diversity_table[2,1],diversity_table[2,2])
	names(allStatTable)[names(allStatTable) == 'tmp'] <- matrixNames[expId]
}

outputStat <- paste(workingPath, "stats-by-plot-454-", otuThr, ".txt", sep = "")
write.table(allStatTable, outputStat, sep="\t", quote=FALSE, row.names=FALSE)


######## non 454 #######	

allStatTable <- NULL
allStatTable$stat <- c("individuals","OTUs","alpha","effective alpha")

for (expId in 1:m) {    	
    communityMatrix <- initNon454ByPlot(expId, otuThr)
    
	source("TurnoverDist.R", local=TRUE)
	# need to put on the top
	source("SampleCounts.R", local=TRUE)

	source("Diversities.R", local=TRUE)

    print(paste(matrixNamesNo454[expId], " (", otuThr, "%) : OTUs = ", ncol(communityMatrix), ", total reads = ", sum(communityMatrix), "; min sample per plot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))
	print(paste("alpha = ", diversity_table[2,1], "; effective alpha = ", diversity_table[2,2], sep=""))
	print(diversity_table)
	
	allStatTable$tmp <- c(sum(communityMatrix),ncol(communityMatrix),diversity_table[2,1],diversity_table[2,2])
	names(allStatTable)[names(allStatTable) == 'tmp'] <- matrixNames[expId]
}

outputStat <- paste(workingPath, "stats-by-plot-traditional-", otuThr, ".txt", sep = "")
write.table(allStatTable, outputStat, sep="\t", quote=FALSE, row.names=FALSE)

