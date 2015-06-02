# maximize diversities with only ? plots
library(vegetarian)

#lev="gamma"
#q=0
#if(!exists("lev")) stop("level of Jost diversity is missing !")
#if(!exists("q")) stop("index of Jost diversity is missing !")


JostDiversity <- function(communityMatrix, level, qN){
	return(d(communityMatrix,lev=level,q=qN))
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
