library(gplots)

if(!exists("sourcePath")) stop("source path to initiate modules is missing !")
if(!exists("workingPath")) stop("working path containing data is missing !")

setwd(sourcePath)

replicate = 1
threshold = 150

# create folder for figures in workingPath
subDir <- "figures"
if (!file.exists(paste(workingPath, subDir, sep=""))) {    
    dir.create(paste(workingPath, subDir, sep=""))    
}

init <- function(expId, otuThr, stringBySubOrPlot) {    
	matrixName <- matrixNames[expId]
	experiment <- experiments[expId]
	
	# e.g. otus97/16S_97_cm.csv
	inputCM <- paste(workingPath, experiments[expId], "/otus", otuThr, "/", experiments[expId], "-", otuThr, ".csv", sep="")
	
	#print(paste("expId = ", expId, ", matrixName = ", matrixName, ", experiment = ", experiment, ", inputCM = ", inputCM, sep=""))  
	
	matrixName <- paste(matrixName, stringBySubOrPlot, sep = "")

	source("InputBySubplot.R", local=TRUE)
	
	source("CommunityTransposeFilter.R", local=TRUE)
		
	return (communityMatrix)
}

# plot_elevations.txt 
#SampleName	Elevation.m		slope.degree	aspect.degree
#Plot1-B	50	18	270
#Plot1-L	50	18	270
#Plot2-C	90	30	60
#Plot2-M	90	30	60
#Plot3-H	160	36	120
#Plot3-N	160	36	120

getElevSampleDist <- function(cmRowNames) { 
    colSample = 1
    colElev = 2
    
    sample_count <- length(cmRowNames)
	print(paste("sample_count = ", sample_count, sep=""))

	inputElevation <- paste(workingPath, "data/plot_elevations.txt", sep="")	
	elev <- read.table(inputElevation, header=TRUE, stringsAsFactors=FALSE)

    # keep elev for easy debug
	elevTmp <- elev
	
   # hard code to fit in given data
    if (sample_count <= 10) {
		# 2 subplots within the plot should be next each other 
		for (i in seq(1,nrow(elevTmp),by=2)) {
			elevTmp[i,colSample] <- unlist(strsplit(as.character(elevTmp[i,colSample]), "-", fixed = TRUE))[1] # hard code for names Plot1-B
			elevTmp[i+1,colSample] <- unlist(strsplit(as.character(elevTmp[i+1,colSample]), "-", fixed = TRUE))[1] 
			if(elevTmp[i,colSample]  != elevTmp[i+1,colSample] ) stop("Two subplots within the same plot should be next each other in the meta data file !") 
		}
	  		
	  	# map subplots meta data to plots
		elevTmp <- elevTmp[-seq(1,nrow(elevTmp), by=2),]		
	} 
	
	# missing data from traditional method
	if (sample_count == 8) {
	  print(paste("Remove sample ", elevTmp[c(7,8),colSample], " from meta data file due to miss data."))
	  elevTmp <- elevTmp[c(-7,-8),]
	} 

    # remove not mapped samples from tmpEnvData
	if ( all(cmRowNames != elevTmp[,colSample]) ) {
		stop(paste("Incorrect plots mapping between CM and meta data :  cm = ", cmRowNames, ", meta data =", elevTmp[,colSample]))
	} 

	elevPlotDist <- dist(elevTmp[,colElev]) # 2nd col Elevation.m
}


######## hard code to deal with traditional methods data ########
# all samples are by plot

# matrixNamesNo454 <-  c("seedlings","trees","invertebrates","birds") 
inputCMNo454 <-  c("seedlings_pilot_species_1_2_0", "", "", "birds_pilot_species_1_2_0")

initNon454ByPlot <- function(expId, otuThr) {    
	experiment <- matrixNamesNo454[expId]
	matrixName <- paste(experiment, sep = "")

    if (experiment == "invertebrates") {
		inputCM <- paste(workingPath, "data/CO1_Invertebrate_Pitfall_1526_OTUs.csv", sep="") # 8 plots
	    cm_subplots <- read.csv(file=inputCM, head=TRUE, sep=",", row.names=1, check.names=FALSE)
        
        print(paste("original data: ", ncol(cm_subplots), " columns = {", paste(colnames(cm_subplots), collapse = ', '), "}, ", nrow(cm_subplots), " rows, sum = ", sum(cm_subplots), ", file = ", inputCM, sep=""))

        communityMatrix <- cm_subplots    
		
		inputCM <- paste(workingPath, "data/CO1_Leaf_Litter_1526_OTUs.csv", sep="") # 10 plots
		cm_subplots <- read.csv(file=inputCM, head=TRUE, sep=",", row.names=1, check.names=FALSE)
        cm_subplots <- cm_subplots[,c(-7,-8)] # remove column 7, 8
        print(paste("after remove column 7, 8: ", ncol(cm_subplots), " columns = {", paste(colnames(cm_subplots), collapse = ', '), "}, ", nrow(cm_subplots), " rows, sum = ", sum(cm_subplots), ", file = ", inputCM, sep=""))  

        communityMatrix <- mergeCommunity(communityMatrix, cm_subplots) 

        print(paste("merged community matrix ", matrixName, ": ", ncol(communityMatrix), " columns = {", paste(colnames(cm_subplots), collapse = ', '), "}, ", nrow(communityMatrix), " rows, sum = ", sum(communityMatrix), sep=""))  

        communityMatrix <- as.matrix(communityMatrix)
		
		
	} else if (experiment == "trees") {
	    inputCM <- paste(workingPath, "data/diameters_pilot_species_1_2_0.csv", sep="")
	    cm_subplots <- read.csv(file=inputCM, head=TRUE, sep=",", row.names=1, check.names=FALSE)
        
        print(paste("original data: ", ncol(cm_subplots), " columns = {", paste(colnames(cm_subplots), collapse = ', '), "}, ", nrow(cm_subplots), " rows, sum = ", sum(cm_subplots), ", file = ", inputCM, sep=""))

        communityMatrix <- cm_subplots    

        inputCM <- paste(workingPath, "data/saplings_pilot_species_1_2_0.csv", sep="")
        cm_subplots <- read.csv(file=inputCM, head=TRUE, sep=",", row.names=1, check.names=FALSE)
        
        print(paste("original data: ", ncol(cm_subplots), " columns = {", paste(colnames(cm_subplots), collapse = ', '), "}, ", nrow(cm_subplots), " rows, sum = ", sum(cm_subplots), ", file = ", inputCM, sep=""))  

        communityMatrix <- mergeCommunity(communityMatrix, cm_subplots) 

        print(paste("merged community matrix ", matrixName, ": ", ncol(communityMatrix), " columns = {", paste(colnames(cm_subplots), collapse = ', '), "}, ", nrow(communityMatrix), " rows, sum = ", sum(communityMatrix), sep=""))  

        communityMatrix <- as.matrix(communityMatrix)
	
	} else {
	  	inputCM <- paste(workingPath, "data/", inputCMNo454[expId], ".csv", sep="")
	    source("InputByPlot.R", local=TRUE)
	}
			
	source("CommunityTransposeFilter.R", local=TRUE)
		
	return (communityMatrix)
}

mergeCommunity <- function(community1, community2) {    
    if(ncol(community1) != ncol(community2)) stop("2 community matrix do not have same columns !")
    if( ! identical(colnames(community1), colnames(community2)) ) stop("2 community matrix do not have same plot !")
    
    comm2 <- merge(community1, community2, by="row.names", all=T)
    comm2[is.na(comm2)] <- 0
    
    # 1st column is "Row.names"
    l <- (length(comm2)-1) / 2
	comm <- data.frame(row.names=comm2[,1])
	for (i in 2:(l+1)) {
		newName <- unlist(strsplit(colnames(comm2)[i], ".", fixed = TRUE))[1]; # hard code for names
		comm[,newName] <- rowSums(comm2[,c(i,i+l)]); # hard code for specific format
	} 
    
    if( (sum(as.matrix(community1)) + sum(as.matrix(community2))) != sum(as.matrix(comm)) ) stop("error during merging !")
    
    return (comm)
}

