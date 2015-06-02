library(xtable)

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

m_comb = 4

getDiversityTable <- function(communityMatrix, matrixName) { 
	source("TurnoverDist.R", local=TRUE)
	# need to put on the top
	source("SampleCounts.R", local=TRUE)

	source("Diversities.R", local=TRUE)

    return(diversity_table)
}

getDiversitiesOfCombPlots <- function(communityMatrix, m_comb) { 
	plotsComb <- combn(rownames(communityMatrix), m_comb)
    d.comb <- matrix(0,nrow=9,ncol=ncol(plotsComb))
    rownames(d.comb) <- paste("$^",qs,"D_\\",levels,"$", sep="") #$^0D_\\gamma$
    
    for (i in 1:ncol(plotsComb)) {
        plots <- plotsComb[,i]
        matrixName <- paste(plots, collapse=" ")
    
        plotId <- which(rownames(communityMatrix) %in% plots)
        
        if (length(plotId) != m_comb) 
			stop(paste("Cannot find plot index from community matrix", plots))
    
		diversity_table <- getDiversityTable(communityMatrix[plotId,], matrixName)

		for (j in  1:9) {
		   d.comb[j,i] <- unlist(diversity_table)[j]
		}
	}
	
	maxD <- matrix(0,nrow=9,ncol=(nrow(communityMatrix)+1))
	rownames(maxD) <- paste("$^",qs,"D_\\",levels,"$", sep="") #$^0D_\\gamma$
	colnames(maxD) <- c(rownames(communityMatrix), "Maximum")
	minD <- matrix(0,nrow=9,ncol=(nrow(communityMatrix)+1))
	rownames(minD) <- paste("$^",qs,"D_\\",levels,"$", sep="") #$^0D_\\gamma$
	colnames(minD) <- c(rownames(communityMatrix), "Minimum")
			
	for (j in  1:9) {
		maxD[j,ncol(maxD)] <- max(d.comb[j,])
		maxId <- which(d.comb[j,] == max(d.comb[j,]))
		plots <- plotsComb[,maxId]
		for (p in plots) {
			plotId <- which(colnames(maxD) %in% p)
			maxD[j,plotId] <- maxD[j,plotId] + 1
		}
		# convert frequency to ratio, last col is max value
		maxD[j,-ncol(maxD)] <- maxD[j,-ncol(maxD)]/max(maxD[j,-ncol(maxD)])  
		
		minD[j,ncol(minD)] <- min(d.comb[j,])
		minId <- which(d.comb[j,] == min(d.comb[j,]))
		plots <- plotsComb[,minId]
		for (p in plots) {
			plotId <- which(colnames(minD) %in% p)
			minD[j,plotId] <- minD[j,plotId] + 1
		}
		# convert frequency to ratio, last col is min value
		minD[j,-ncol(minD)] <- minD[j,-ncol(minD)]/max(minD[j,-ncol(minD)]) 
		
		print(paste("max", rownames(d.comb)[j], "index =", maxId, "; min index =", minId))
	}

#	maxD[,ncol(maxD)] <- formatC(signif(maxD[,ncol(maxD)],digits=5), digits=5,format="fg", flag="#")
#	minD[,ncol(minD)] <- formatC(signif(minD[,ncol(minD)],digits=5), digits=5,format="fg", flag="#")
#	maxD[maxD==0] <- ""
#	minD[minD==0] <- ""

	colnames(maxD) <- gsub("CM30C30", "Plot9", colnames(maxD), ignore.case = T)
	colnames(maxD) <- gsub("LB1", "Plot10", colnames(maxD), ignore.case = T)
	colnames(minD) <- gsub("CM30C30", "Plot9", colnames(minD), ignore.case = T)
	colnames(minD) <- gsub("LB1", "Plot10", colnames(minD), ignore.case = T)

	print(maxD)
	print(minD)

	list("maxD" = maxD, "minD" = minD)
}

######## 454 #######	

for (expId in 1:n) {	
    communityMatrix <- init(expId, otuThr, "-by-plot")
    
    print(paste(matrixNames[expId], " (", otuThr, "%) : OTUs = ", ncol(communityMatrix), ", total reads = ", sum(communityMatrix), "; min sample per plot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))
		    
	######## combinations of 4 #######
	rangeD <- getDiversitiesOfCombPlots(communityMatrix, m_comb) 
	
	maxD <- rangeD$maxD
	
#	minD <- rangeD$minD
#	print(xtable(maxD),sanitize.text.function=function(x){x})
#	print(xtable(minD),sanitize.text.function=function(x){x})

	# removed max value
	maxD <- maxD[,-ncol(maxD)]
		
	if (expId==1) {
		gamma_0 <- matrix(0,nrow=nrow(communityMatrix),ncol=(n+m))
		rownames(gamma_0) <- rownames(communityMatrix)
		colnames(gamma_0) <- c(matrixNames, matrixNamesNo454)

		beta_1 <- matrix(0,nrow=nrow(communityMatrix),ncol=(n+m))
		rownames(beta_1) <- rownames(communityMatrix)
		colnames(beta_1) <- c(matrixNames, matrixNamesNo454)
		
		rownames(gamma_0) <- gsub("CM30C30", "Plot9", rownames(gamma_0), ignore.case = T)
		rownames(gamma_0) <- gsub("LB1", "Plot10", rownames(gamma_0), ignore.case = T)
		rownames(beta_1) <- gsub("CM30C30", "Plot9", rownames(beta_1), ignore.case = T)
		rownames(beta_1) <- gsub("LB1", "Plot10", rownames(beta_1), ignore.case = T)	
	}

	if(all(tolower(rownames(gamma_0)) != tolower(colnames(maxD)))) 
		stop(paste("Plot name does not match in ", matrixNames[expId]))

	# gamma_0	
	gamma_0[,expId] <- maxD[1,]
	# beta_1	
	beta_1[,expId] <- maxD[6,]
}

######## non 454 #######	

for (expId in 1:m) {    	
    communityMatrix <- initNon454ByPlot(expId, otuThr)
        
    print(paste(matrixNamesNo454[expId], " (", otuThr, "%) : OTUs = ", ncol(communityMatrix), ", total reads = ", sum(communityMatrix), "; min sample per plot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))
		    
	######## combinations of 4 #######
	rangeD <- getDiversitiesOfCombPlots(communityMatrix, m_comb) 

	maxD <- rangeD$maxD
#	minD <- rangeD$minD
#	print(xtable(maxD),sanitize.text.function=function(x){x})
#	print(xtable(minD),sanitize.text.function=function(x){x})

	# removed max value
	maxD <- maxD[,-ncol(maxD)]
		
	if (expId==3) { # hard code for missing data
		maxD <- cbind(maxD, Plot7=rep(0,nrow(maxD)))
		maxD <- cbind(maxD, Plot8=rep(0,nrow(maxD)))
		maxD <- maxD[,c(1:6,9,10,7:8)]
	} 
	
	if(all(tolower(rownames(gamma_0)) != tolower(colnames(maxD)))) 
		stop(paste("Plot name does not match in ", matrixNamesNo454[expId]))

	# gamma_0	
	gamma_0[,(expId+n)] <- maxD[1,]
	# beta_1	
	beta_1[,(expId+n)] <- maxD[6,]
}

######## file #######
outputRank <- paste(workingPath, "prob-4plots-table-gamma0.txt", sep = "")
write.table(gamma_0, outputRank, sep="\t", row.names=T)

outputRank <- paste(workingPath, "prob-4plots-table-beta1.txt", sep = "")
write.table(beta_1, outputRank, sep="\t", row.names=T)

print(xtable(gamma_0, caption = "The probability of having maximum Jost diversities of all possible combinations of 4 plots", 
	label = "tab:4plotsMaxGamma0"), sanitize.text.function = function(x){x})
	
print(xtable(beta_1, caption = "The probability of having maximum Jost diversities of all possible combinations of 4 plots", 
	label = "tab:4plotsMaxBeta1"), sanitize.text.function = function(x){x})
	
	
