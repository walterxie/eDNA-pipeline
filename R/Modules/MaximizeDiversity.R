# maximize diversities with only ? plots
library(vegetarian)

JostDiversity <- function(communityMatrix, level, qN){
	return(d(communityMatrix,lev=level,q=qN))
}

plotHeatmapMaxRemDivRank <- function(rank.elv, fname, title) {
  breaks.rank <- round(seq(1, nrow(rank.elv), length.out = 5), digits = 0)
  
  rank.melt <- melt(rank.elv, id=c("Samples","Elevation","ForestType"))
  rank.melt$Samples <- factor(rank.melt$Samples, levels=unique(rank.melt$Samples))
  
  p <- ggplot(rank.melt, aes(x=variable, y=Samples)) + 
    geom_tile(aes(fill=value)) + 
    scale_fill_gradient(na.value="transparent", low="white", high="steelblue", name="rank", breaks=breaks.rank) +
    ylab("Plot (sorted by elevation)") + ggtitle(title) +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
  
  pdf(file.path(workingPath, figDir, fname), width=6, height=6)	
  print(p)
  invisible(dev.off()) 
}



# Note: communityMatrix must be by plots
# return 3 columns: rank, diversity, removedSites
#e.g.
#rank	diversity	removedSites
#28	1845.785714	CM30c39
#27	1875.888889	CM30c44
#26	1899.653846	CM31a5
#25	1925.08	CM30c46
#24	1950.708333	Plot1
# greedy algorithm, always choose the 1st plot for multi-results 
getMaxRemainedDiversity <- function(communityMatrix, level, qN){
	tmpCM <- communityMatrix
	removedSites <- c()

	jostDiversity <- JostDiversity(communityMatrix, lev, q)
	# add orignal diversity
	maxRemainedDiversities <- c(jostDiversity)

	for (ra in nrow(communityMatrix):2) {
		if(ra != nrow(tmpCM)) stop(paste("incorrect nrow(tmpCM) ", nrow(tmpCM), " !=  ra ", ra, sep=""))
	
		maxRemainedDiv = 0
		removedSiteId = 0

		for (siteId in 1:nrow(tmpCM)) {
			tmpDiv <- JostDiversity(tmpCM[-siteId,], lev, q)
			if (tmpDiv > maxRemainedDiv) {
				maxRemainedDiv = tmpDiv
				removedSiteId = siteId
			}
		
#			print(paste("siteId = ", siteId, ", remove siteId = ", removedSiteId, ", maxRemainedDiv = ", maxRemainedDiv, sep=""))  
		}

		if(maxRemainedDiv < 1 | removedSiteId < 1) stop(paste("incorrect max remained ", lev, "(", q, ") = ", maxRemainedDiv, ", when sites = ", nrow(tmpCM), sep=""))
	
		print(paste("remove site ", rownames(tmpCM)[removedSiteId], ", max remained ", lev, "(", q, ") = ", maxRemainedDiv, sep=""))  
		removedSites <- c(removedSites, rownames(tmpCM)[removedSiteId])
		maxRemainedDiversities <- c(maxRemainedDiversities, maxRemainedDiv)
	
		tmpCM <- tmpCM[-removedSiteId,]       
	}
	removedSites <- c(removedSites, rownames(tmpCM)[1]) # the last

	if(length(removedSites) != length(maxRemainedDiversities)) 
		stop(paste("length of maxRemainedDiversities ", length(maxRemainedDiversities), " - 1 != sites ", length(maxDiv$rank), sep=""))

	# 1st removed (biggest rank) maxDiv$rank is the least important
	maxDiv <- data.frame(row.names=removedSites, rank=nrow(communityMatrix):1, diversity=maxRemainedDiversities) 
	maxDiv <- maxDiv[order(rownames(maxDiv)),]

	return(maxDiv)
}

# give same rank to the same remaining
getMaxRemainedDiversity2 <- function(communityMatrix, level, qN){
	tmpCM <- communityMatrix[order(rownames(communityMatrix)),]

	ranks<-c()
	removedSites <- c()
	maxRemainedDiversities<-c()
	while(nrow(tmpCM) > 1) {		
		divs <- c()
		for (siteId in 1:nrow(tmpCM)) {
			tmpDiv <- JostDiversity(tmpCM[-siteId,], lev, q)
			divs <- c(divs, tmpDiv)
		
			print(paste("siteId = ", siteId, ", tmpDiv = ", tmpDiv, sep="")) 
		}
		maxRemainedDiv = max(divs)	
		removedSiteId = which(divs==maxRemainedDiv)

		if(maxRemainedDiv < 0 | length(removedSiteId) < 1) 
			stop(paste("incorrect max remained ", lev, "(", q, ") = ", maxRemainedDiv, ", when sites = ", nrow(tmpCM), sep=""))
	
		print(paste("remove site ", rownames(tmpCM)[removedSiteId], ", max remained ", lev, "(", q, ") = ", maxRemainedDiv, sep=""))  
		
		ranks<-c(ranks, rep(nrow(tmpCM), length(removedSiteId)))
		removedSites <- c(removedSites, rownames(tmpCM)[removedSiteId])
		maxRemainedDiversities <- c(maxRemainedDiversities, rep(maxRemainedDiv, length(removedSiteId)))
	
		tmpCM <- tmpCM[-removedSiteId,]       
	}
	
	if (nrow(tmpCM) == 1) {
		ranks<-c(ranks, 1)
		removedSites <- c(removedSites, rownames(tmpCM)[1]) # the last
		maxRemainedDiversities <- c(maxRemainedDiversities, 0)
	}

	if(length(ranks) != length(maxRemainedDiversities)) 
		stop(paste("length of maxRemainedDiversities ", length(maxRemainedDiversities), " - 1 != sites ", length(ranks), sep=""))

	# 1st removed (biggest rank) maxDiv$rank is the least important
	maxDiv <- data.frame(row.names=removedSites, rank=ranks, diversity=maxRemainedDiversities)
	maxDiv <- maxDiv[order(rownames(maxDiv)),]

	return(maxDiv)
}

# get diversity_table from Diversities.R
getDiversityTable <- function(communityMatrix, matrixName) { 
	source("Modules/TurnoverDist.R", local=TRUE)
	# need to put on the top
	source("Modules/SampleCounts.R", local=TRUE)

	source("Modules/Diversities.R", local=TRUE)

    return(diversity_table)
}

# maximum and minimum Jost diversity of all possible combinations of m_comb plots
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
