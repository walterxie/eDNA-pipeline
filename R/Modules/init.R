
if(!exists("sourcePath")) stop("source path to initiate modules is missing !")
if(!exists("workingPath")) stop("working path containing data is missing !")

#setwd(sourcePath)
#verbose <- FALSE

replicate = 1
threshold = 150

if (verbose) 
cat(paste("\nConfig : data has", replicate, "replicate per sample, report", threshold, "most abundant OTUs.\n"))

# create folder for figures in workingPath
#figDir <- "figures"
#if (!file.exists(paste(workingPath, figDir, sep=""))) {    
#    dir.create(paste(workingPath, figDir, sep=""))    
#}

######## load 454 community matrix #######
# expId = 1:6
init <- function(expId, otuThr, stringBySubOrPlot) {    
	matrixName <- matrixNames[expId]
	
	# e.g. otus97/16S_97_cm.csv
	inputCM <- paste(workingPath, matrixName, "/otus", otuThr, "/", matrixName, "-", otuThr, ".csv", sep="")
	
	#print(paste("expId = ", expId, ", matrixName = ", matrixName, ", inputCM = ", inputCM, sep=""))  
	
	matrixName <- paste(matrixName, stringBySubOrPlot, sep = "")

	source("Modules/InputBySubplot.R", local=TRUE)
	
	source("Modules/CommunityTransposeFilter.R", local=TRUE)
	
	if (rmSingleton) {
		singletons <- which(colSums(communityMatrix)==1)
		print(paste("Remove", length(singletons) ,"singletons from ", matrixName, " ! "))
		communityMatrix <- communityMatrix[,-singletons]
		rm(singletons)		
	}
	
	return (communityMatrix)
}


######## get min size all sites for rdiversityTable each gene #######
getMinSizeAllSites <- function(expId, otuThre) {
    minSample = 99999

	communityMatrix <- init(expId, otuThre, stringBySubOrPlot)
	minSample <- min(minSample, min(rowSums(communityMatrix)))  		

	sampleSize <- minSample#floor(minSample/100) * 100 # make sure subsampling run  
	print(paste("otuThre =", otuThre, ", min sample size per site =", minSample, ", take sample size per site =", sampleSize)) 
    
    return (sampleSize)
}

######## beta 1 -1 #######
beta1 <- function(communityMatrix) {    
	# including diagonal
    rowsNum <- nrow(communityMatrix) * (nrow(communityMatrix) + 1) / 2  		
	d.beta1 <- matrix(0,nrow=nrow(communityMatrix),ncol=nrow(communityMatrix))
	colnames(d.beta1) <- c(rownames(communityMatrix))
    rownames(d.beta1) <- c(rownames(communityMatrix))
	count=0
	for(i in 1:nrow(communityMatrix)){
		for(k in i:nrow(communityMatrix)){
			count=count+1				
			d.beta1[k,i] <- d(communityMatrix[c(i,k),],lev="beta",q=1)-1				
		}
	}
	if (count != rowsNum) stop("incorrect pairwise comparisons !")
	
	return (as.dist(d.beta1))
}
######## alpha1 #######
alpha1 <- function(communityMatrix) {    
	# including diagonal
    m.alpha1 <- matrix(0,nrow=nrow(communityMatrix),ncol=1)	
    rownames(m.alpha1) <- c(rownames(communityMatrix))
	for(i in 1:nrow(communityMatrix)){				
		m.alpha1[i,1] <- d(communityMatrix[i,],lev="gamma",q=1)				
	}
	
	return (m.alpha1) # 1 col matrix
}
######## get table in the format of "corr (sign)" #######
# Input: corr.sign.matrix is a matrix having same row and col names, 
# lower triangle is correlations (or equivalent), upper triangle is significance
# Output corr.sign.table is in the format of "corr (sign)"
getCorrSignTable <- function(corr.sign.matrix, digits=3) {
	m.corr <- corr.sign.matrix
	m.corr[upper.tri(m.corr)] <- 0
	m.corr <- formatC(signif(m.corr,digits=digits), digits=digits,format="fg", flag="#")
	m.sign <- t(corr.sign.matrix)
	m.sign[upper.tri(m.sign)] <- 0
	m.sign <- formatC(signif(m.sign,digits=digits), digits=digits,format="fg", flag="#")

	corr.sign.table <- matrix( paste(m.corr, " (", m.sign, ")", sep=""), nrow=nrow(m.corr), dimnames=dimnames(m.corr) )
	corr.sign.table[corr.sign.table=="0 (0)"] <- ""

	corr.sign.table <- corr.sign.table[-1,-ncol(corr.sign.table)]
}


######## elevations #######
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
######## elevations #######
getElevSample <- function(cmRowNames) { 
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

    elevPlot <- data.frame(row.names=cmRowNames)
    elevPlot$Elevation.m <- elevTmp[,colElev]
	return(elevPlot) # 2nd col Elevation.m
}

######## hard code to deal with traditional methods data ########
# all samples are by plot
# matrixNamesNo454 <-  c("seedlings","trees","invertebrates","birds") 
inputCMNo454 <-  c("seedlings_pilot_species_1_2_0", "", "", "birds_pilot_species_1_2_0")
# expId = 1:4
initNon454ByPlot <- function(expId, otuThr) {    
	matrixName <- matrixNamesNo454[expId]
	
    if (matrixName == "invertebrates") {
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
		
		
	} else if (matrixName == "trees") {
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
	    source("Modules/InputByPlot.R", local=TRUE)
	}
			
	source("Modules/CommunityTransposeFilter.R", local=TRUE)
	
	return (communityMatrix)
}
######## merge 2 community matrix #######
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

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
