# "-by-plot" trigger merge 2 subplots columns
library(vegan)
library(vegetarian)

beta1 <- function(communityMatrix) {    
	####### beta 1 ########
	# including diagonal
    rowsNum <- nrow(communityMatrix) * (nrow(communityMatrix) + 1) / 2  		
	d.beta1 <- matrix(0,nrow=nrow(communityMatrix),ncol=nrow(communityMatrix))
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

# change config below
sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
setwd(sourcePath)

workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
experiments <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name
matrixNames <-  experiments

matrixNamesNo454 <-  c("seedlings","trees","inverts","birds") # need expId correct for "birds","seedlings" 

n <- length(matrixNames) 
m <- length(matrixNamesNo454)

otuThr = 97

source("init.R", local=TRUE)

######## 454 #######
beta1DisList <- list()

for (expId in 1:n) {	
    communityMatrix <- init(expId, otuThr, "-by-plot")

	####### beta 1 ########
	# (45 pairs for 10 plots)
    d.beta1 <- beta1(communityMatrix)
    
    beta1DisList[[ expId ]] <- d.beta1
	    
    print(paste(matrixNames[expId], " (", otuThr, "%) : total = ", sum(communityMatrix), "; min sample per subplot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))
	
}

######## non 454 #######	
beta1No454DisList <- list()

for (expId in 1:m) {    	
    communityMatrix <- initNon454ByPlot(expId, otuThr)
    
	####### beta 1 ########
	# (45 pairs for 10 plots)
    d.beta1 <- beta1(communityMatrix)
    
    beta1No454DisList[[ expId ]] <- d.beta1
       
    print(paste(matrixNamesNo454[expId], " (", otuThr, "%) : total = ", sum(communityMatrix), "; min sample per subplot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))

}

####### mantel correlations (and p-values) between the different communities #######

m.mantel <- matrix(0,nrow=(n+m),ncol=(n+m))
m.signif <- matrix(0,nrow=(n+m),ncol=(n+m))
colnames(m.mantel) <- c(matrixNames, matrixNamesNo454)
rownames(m.mantel) <- c(matrixNames, matrixNamesNo454)

euc.dist <- matrix(0,nrow=(n+m),ncol=(n+m))
colnames(euc.dist) <- c(matrixNames, matrixNamesNo454)
rownames(euc.dist) <- c(matrixNames, matrixNamesNo454)

for (expId1 in 1:(n+m-1)) {
	for (expId2 in (expId1+1):(n+m)) {	     
	     if (expId1 <= n) {
	         d.comm1 <- beta1DisList[[ expId1 ]]
	         exp.comm1 <- matrixNames[expId1]
	     } else {
	         d.comm1 <- beta1No454DisList[[ expId1-n ]]
	         exp.comm1 <- matrixNamesNo454[expId1-n]
	     }
	     
	     if (expId2 <= n) {
	         d.comm2 <- beta1DisList[[ expId2 ]]
	         exp.comm2 <- matrixNames[expId2]
	     } else {
	         d.comm2 <- beta1No454DisList[[ expId2-n ]]
	         exp.comm2 <- matrixNamesNo454[expId2-n]
	     }
	     
	     
	     if (length(d.comm1) > length(d.comm2)) {
			m.comm <- as.matrix(d.comm1) 
			d.comm1 <- as.dist(m.comm[c(-7,-8),c(-7,-8)]) 
	     } else if (length(d.comm1) < length(d.comm2)) {
			m.comm <- as.matrix(d.comm2) 
			d.comm2 <- as.dist(m.comm[c(-7,-8),c(-7,-8)]) 	     
	     }
	     
	     mantel.comm <- mantel(d.comm1, d.comm2, permutations=4999) 
	     
	     m.mantel[expId2, expId1] <- mantel.comm$statistic
	     m.mantel[expId1, expId2] <- mantel.comm$statistic
	     m.signif[expId2, expId1] <- mantel.comm$signif
	     m.signif[expId1, expId2] <- mantel.comm$signif
	     
	     print(paste(exp.comm1, " vs. ", exp.comm2, ", mantel statistic r = ", round(mantel.comm$statistic, 3), "; significance = ", mantel.comm$signif, "; permutations = ", mantel.comm$permutations, sep=""))
	     
	     #euc.dist.comm <- dist(rbind(as.vector(d.comm1), as.vector(d.comm2)), method = "manhattan")
	     #euc.dist[expId2, expId1] <- euc.dist.comm
	     #euc.dist[expId1, expId2] <- euc.dist.comm
	}
}

for (expId1 in 1:(n+m)) {
    m.mantel[expId1, expId1] <- 1
}

###### print table #####
df.mantel <- data.frame(m.mantel, check.names=FALSE)
for (expId1 in 1:(n+m-1)) {
	for (expId2 in (expId1+1):(n+m)) {	     
        df.mantel[expId1, expId2] <- c("")
        tmp <- paste(round(m.mantel[expId2, expId1], 3), " (", m.signif[expId2, expId1], ")", sep="")
        df.mantel[expId2, expId1] <- c(tmp)
	}
}
for (expId1 in 1:(n+m)) {
    df.mantel[expId1, expId1] <- c("")
}
df.mantel <- df.mantel[-1,-10]

mantel.xtable <- xtable(df.mantel)
print(mantel.xtable, sanitize.text.function = function(x){x})


####### graph #####
# Classical MDS
d.mantel <- 1-m.mantel # euclidean distances between the rows
fit <- cmdscale(d.mantel,eig=TRUE, k=2) # k is the number of dim

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]

pdf(paste(workingPath, "figures/mds-pairewise-cm-", otuThr, ".pdf", sep = ""), width=5, height=5)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS",	 type="n")
text(x, y, labels = row.names(m.mantel), cex=.7)
invisible(dev.off()) 

m.mantel1 <- m.mantel[-10,-10] # no birds
d.mantel <- 1-m.mantel1 # euclidean distances between the rows
fit <- cmdscale(d.mantel,eig=TRUE, k=2) # k is the number of dim

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]

pdf(paste(workingPath, "figures/mds-no-birds-pairewise-cm-", otuThr, ".pdf", sep = ""), width=5, height=5)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS",	 type="n")
text(x, y, labels = row.names(m.mantel1), cex=.7)
invisible(dev.off()) 