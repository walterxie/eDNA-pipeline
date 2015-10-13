

# change config below
#sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
#setwd(sourcePath)
#workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
#matrixNames <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name   
#levels = rep(c("gamma","alpha","beta"),3)
#qs = rep(0:2,each=3)

if(!exists("tableFile")) stop("table file is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("isPlot")) isPlot = TRUE # by plot
if(!exists("otuThr")) otuThr = 97
if(!exists("levels")) levels = rep(c("gamma","alpha","beta"),3)
if(!exists("qs")) qs = rep(0:2,each=3)

n <- length(matrixNames) 

source("Modules/init.R", local=TRUE)

divRowNames <- c("$^0D_\\gamma$","$^0D_\\alpha$","$^0D_\\beta$","$^1D_\\gamma$","$^1D_\\alpha$","$^1D_\\beta$",
				"$^2D_\\gamma$","$^2D_\\alpha$","$^2D_\\beta$")
otusRowNames <- c("Reads", "OTUs", "Singleton", "Coupleton") # cannot get unique reads from CM

## gene by gene, all plots 
allDiv <- matrix(0,nrow=length(divRowNames),ncol=n)
colnames(allDiv) <- matrixNames
rownames(allDiv) <- divRowNames

# OTU stats with total
allOTUs <- matrix(0,nrow=length(otusRowNames),ncol=(n+1))
colnames(allOTUs) <- c(matrixNames, "Total")
rownames(allOTUs) <- otusRowNames

######## summary #######

cat("Table: data summary: rmSingleton =", rmSingleton, ", isPlot =", isPlot, ", otuThr =", otuThr, "/n") 

for (expId in 1:n) {
	# "-by-plot" trigger merge 2 subplots columns	
    communityMatrix <- init(expId, otuThr, "-by-plot")
    # for Diversities module
    matrixName <- paste(matrixNames[expId], "-by-plot", sep = "")
    
	source("Modules/TurnoverDist.R", local=TRUE)
	# need to put on the top
	source("Modules/SampleCounts.R", local=TRUE)

	source("Modules/Diversities.R", local=TRUE)

	print(paste(matrixNames[expId], " (", otuThr, "%) : OTUs = ", ncol(communityMatrix), ", total reads = ", sum(communityMatrix), "; min sample per plot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))
	print(paste("alpha = ", diversity_table[2,1], "; effective alpha = ", diversity_table[2,2], sep=""))
#	print(diversity_table)
    
    for (i in 1:length(unlist(diversity_table))) {
		allDiv[i,expId] <- unlist(diversity_table)[i]
    }
    allOTUs[1,expId] <- sum(communityMatrix)
    allOTUs[2,expId] <- ncol(communityMatrix)
    allOTUs[3,expId] <- sum(colSums(communityMatrix)==1)
    allOTUs[4,expId] <- sum(colSums(communityMatrix)==2)

	## plot by plot
	if (expId==1) {	
		# gamme0
		allReadsPerPlot <- matrix(0,nrow=nrow(communityMatrix),ncol=(n+1))
		colnames(allReadsPerPlot) <- c(matrixNames, "Total")
		rownames(allReadsPerPlot) <- rownames(communityMatrix)	
		# OTUs per plot
		allOTUsPerPlot <- matrix(0,nrow=nrow(communityMatrix),ncol=(n+1))
		colnames(allOTUsPerPlot) <- c(matrixNames, "Total")
		rownames(allOTUsPerPlot) <- rownames(communityMatrix)
		# Shannon
		allGamma1PerPlot <- matrix(0,nrow=nrow(communityMatrix),ncol=n)
		colnames(allGamma1PerPlot) <- matrixNames
		rownames(allGamma1PerPlot) <- rownames(communityMatrix)
	}

	allReadsPerPlot[,expId] <- rowSums(communityMatrix)
	for (i in 1:nrow(communityMatrix)) {		
		allOTUsPerPlot[i,expId] <- sum(communityMatrix[i,]>0)
		allGamma1PerPlot[i,expId] <- d(communityMatrix[i,],lev="gamma",q=1)
	}
}
## gene by gene, all plots
for (i in 1:ncol(allDiv)) {
	allDiv[,i] <- as.numeric(allDiv[,i])
} 
allDiv <- round(allDiv, 2)
allDiv <- format(allDiv, big.mark=",", scientific=F)

print(xtable(allDiv, caption = "Table of Jost diversities for eDNA data sets", label = "tab:biodiveDNA", caption.placement = "top"), 
		sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

allOTUs[,(n+1)] <- rowSums(allOTUs)
allOTUs <- format(allOTUs, big.mark=",", scientific=F)

print(xtable(allOTUs, caption = "Table of OTUs for eDNA data sets", label = "tab:biodiveDNA", caption.placement = "top"), 
		sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

## plot by plot
allReadsPerPlot[,(n+1)] <- rowSums(allReadsPerPlot)
allReadsPerPlot <- format(allReadsPerPlot, big.mark=",", scientific=F)

print(xtable(allReadsPerPlot, caption = "Table of sequences per plot for eDNA data sets", label = "tab:readseDNAPerPlot", caption.placement = "top"), 
		sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

allOTUsPerPlot[,(n+1)] <- rowSums(allOTUsPerPlot)
allOTUsPerPlot <- format(allOTUsPerPlot, big.mark=",", scientific=F)

print(xtable(allOTUsPerPlot, caption = "Table of OTUs per plot for eDNA data sets", label = "tab:otuseDNAPerPlot", caption.placement = "top"), 
		sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

for (i in 1:ncol(allGamma1PerPlot)) {
	allGamma1PerPlot[,i] <- as.numeric(allGamma1PerPlot[,i])
}
allGamma1PerPlot <- round(allGamma1PerPlot, 2)
allGamma1PerPlot <- format(allGamma1PerPlot, big.mark=",", scientific=F)

print(xtable(allGamma1PerPlot, caption = "Table of Shannon index per plot for eDNA data sets", label = "tab:gamma1eDNAPerPlot", caption.placement = "top"), 
		sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)


######## non 454 #######	

## gene by gene, all plots 
allDiv <- matrix(0,nrow=length(divRowNames),ncol=m)
colnames(allDiv) <- matrixNamesNo454
rownames(allDiv) <- divRowNames

# OTU stats with total
allOTUs <- matrix(0,nrow=length(otusRowNames),ncol=(m+1))
colnames(allOTUs) <- c(matrixNamesNo454, "Total")
rownames(allOTUs) <- otusRowNames

for (expId in 1:m) {    	
    communityMatrix <- initNon454ByPlot(expId, otuThr)
    matrixName <- matrixNamesNo454[expId]
    
	source("Modules/TurnoverDist.R", local=TRUE)
	# need to put on the top
	source("Modules/SampleCounts.R", local=TRUE)

	source("Modules/Diversities.R", local=TRUE)

	print(paste(matrixNames[expId], " (", otuThr, "%) : OTUs = ", ncol(communityMatrix), ", total reads = ", sum(communityMatrix), "; min sample per plot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))
	print(paste("alpha = ", diversity_table[2,1], "; effective alpha = ", diversity_table[2,2], sep=""))
#	print(diversity_table)
    
    for (i in 1:length(unlist(diversity_table))) {
		allDiv[i,expId] <- unlist(diversity_table)[i]
    }
    allOTUs[1,expId] <- sum(communityMatrix)
    allOTUs[2,expId] <- ncol(communityMatrix)
    allOTUs[3,expId] <- sum(colSums(communityMatrix)==1)
    allOTUs[4,expId] <- sum(colSums(communityMatrix)==2)

	## plot by plot
	if (expId==1) {	
		# gamme0
		allReadsPerPlot <- matrix(0,nrow=nrow(communityMatrix),ncol=(m+1))
		colnames(allReadsPerPlot) <- c(matrixNamesNo454, "Total")
		rownames(allReadsPerPlot) <- rownames(communityMatrix)	
		# OTUs per plot
		allOTUsPerPlot <- matrix(0,nrow=nrow(communityMatrix),ncol=(m+1))
		colnames(allOTUsPerPlot) <- c(matrixNamesNo454, "Total")
		rownames(allOTUsPerPlot) <- rownames(communityMatrix)
		# Shannon
		allGamma1PerPlot <- matrix(0,nrow=nrow(communityMatrix),ncol=m)
		colnames(allGamma1PerPlot) <- matrixNamesNo454
		rownames(allGamma1PerPlot) <- rownames(communityMatrix)
	}

	if (expId==3) {
		# hard code for missing data
		allReadsPerPlot[c(1:6,9,10),expId] <- rowSums(communityMatrix)
		for (i in 1:6) {		
			allOTUsPerPlot[i,expId] <- sum(communityMatrix[i,]>0)
			allGamma1PerPlot[i,expId] <- d(communityMatrix[i,],lev="gamma",q=1)
		}
		for (i in 7:8) {			
			allOTUsPerPlot[(i+2),expId] <- sum(communityMatrix[i,]>0)
			allGamma1PerPlot[(i+2),expId] <- d(communityMatrix[i,],lev="gamma",q=1)
		}
		
	} else {
		allReadsPerPlot[,expId] <- rowSums(communityMatrix)
		for (i in 1:nrow(communityMatrix)) {		
			allOTUsPerPlot[i,expId] <- sum(communityMatrix[i,]>0)
			allGamma1PerPlot[i,expId] <- d(communityMatrix[i,],lev="gamma",q=1)
		}
	}
}
## gene by gene, all plots
for (i in 1:ncol(allDiv)) {
	allDiv[,i] <- as.numeric(allDiv[,i])
} 
allDiv <- round(allDiv, 2)
allDiv <- format(allDiv, big.mark=",", scientific=F)

print(xtable(allDiv, caption = "Table of Jost diversities for traditional data sets", label = "tab:biodivTrad", caption.placement = "top"), 
		sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

allOTUs[,(m+1)] <- rowSums(allOTUs)
allOTUs <- format(allOTUs, big.mark=",", scientific=F)

print(xtable(allOTUs, caption = "Table of OTUs for traditional data sets", label = "tab:biodivTrad", caption.placement = "top"), 
		sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

## plot by plot
allReadsPerPlot[,(m+1)] <- rowSums(allReadsPerPlot)
allReadsPerPlot <- format(allReadsPerPlot, big.mark=",", scientific=F)

print(xtable(allReadsPerPlot, caption = "Table of sequences per plot for traditional data sets", label = "tab:readsTradPerPlot", caption.placement = "top"), 
		sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

allOTUsPerPlot[,(m+1)] <- rowSums(allOTUsPerPlot)
allOTUsPerPlot <- format(allOTUsPerPlot, big.mark=",", scientific=F)

print(xtable(allOTUsPerPlot, caption = "Table of OTUs per plot for traditional data sets", label = "tab:otusTradPerPlot", caption.placement = "top"), 
		sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

for (i in 1:ncol(allGamma1PerPlot)) {
	allGamma1PerPlot[,i] <- as.numeric(allGamma1PerPlot[,i])
}
allGamma1PerPlot <- round(allGamma1PerPlot, 2)
allGamma1PerPlot <- format(allGamma1PerPlot, big.mark=",", scientific=F)

print(xtable(allGamma1PerPlot, caption = "Table of Shannon index per plot for traditional data sets", label = "tab:gamma1TradPerPlot", caption.placement = "top"), 
		sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

