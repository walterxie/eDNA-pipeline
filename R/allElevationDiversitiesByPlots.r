# "-by-plot" trigger merge 2 subplots columns
library(RColorBrewer)
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

matrixNamesNo454 <-  c("seedlings","trees","invertebrates","birds") # need expId correct for "birds","seedlings" 

n <- length(matrixNames) 
mypalette <- c("red", "orange", "green", "purple", "blue", "brown")
myshape <- seq(0, (0 + n-1))

m <- length(matrixNamesNo454)
mypalette2 <- brewer.pal(m,"Dark2")
myshape2 <- seq(0, (0 + m-1))

otuThr = 97

source("init.R", local=TRUE)

######## 454 #######
whittakerDisList <- list()
beta1DisList <- list()
hornMorisitaDisList <- list()
elevPlotDistList <- list()

for (expId in 1:n) {	
    communityMatrix <- init(expId, otuThr, "-by-plot")
    elevPlotDist <- getElevSampleDist(rownames(communityMatrix)) 
    
	######## Whittaker's beta #######	
	d.brayBin <- vegdist(communityMatrix, method="bray", binary=TRUE)
	#mantel.brayBin <- mantel(elevPlotDist, d.brayBin, permutations=4999)
    
    whittakerDisList[[ expId ]] <- d.brayBin
	elevPlotDistList[[ expId ]] <- elevPlotDist	
	
	####### beta 1 ########
	# (45 pairs for 10 plots)
    d.beta1 <- beta1(communityMatrix)
    #mantel.beta1 <- mantel(elevPlotDist, d.beta1, permutations=4999)
    
    beta1DisList[[ expId ]] <- d.beta1
	
	####### Horn-Morisita ########
    d.hornMorisita <- vegdist(communityMatrix, method="horn", binary=FALSE)
	#mantel.hornMorisita <- mantel(elevPlotDist, d.hornMorisita, permutations=4999)
    
    hornMorisitaDisList[[ expId ]] <- d.hornMorisita
    
    print(paste(matrixNames[expId], " (", otuThr, "%) : total = ", sum(communityMatrix), "; min sample per subplot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))
	
}

######## non 454 #######	
whittakerNo454DisList <- list()
beta1No454DisList <- list()
hornMorisitaNo454DisList <- list()
elevPlotDistNo454List <- list()

for (expId in 1:m) {    	
    communityMatrix <- initNon454ByPlot(expId, otuThr)
    elevPlotDist <- getElevSampleDist(rownames(communityMatrix)) 
    
	######## Whittaker's beta #######	
	d.brayBin <- vegdist(communityMatrix, method="bray", binary=TRUE)
	#mantel.brayBin <- mantel(elevPlotDist, d.brayBin, permutations=4999)
    
    whittakerNo454DisList[[ expId ]] <- d.brayBin
	elevPlotDistNo454List[[ expId ]] <- elevPlotDist
	
	####### beta 1 ########
	# (45 pairs for 10 plots)
    d.beta1 <- beta1(communityMatrix)
    #mantel.beta1 <- mantel(elevPlotDist, d.beta1, permutations=4999)
    
    beta1No454DisList[[ expId ]] <- d.beta1
    
    ####### Horn-Morisita ########
    d.hornMorisita <- vegdist(communityMatrix, method="horn", binary=FALSE)
	#mantel.hornMorisita <- mantel(elevPlotDist, d.hornMorisita, permutations=4999)
    
    hornMorisitaNo454DisList[[ expId ]] <- d.hornMorisita
    
    print(paste(matrixNamesNo454[expId], " (", otuThr, "%) : total = ", sum(communityMatrix), "; min sample per subplot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))

}

####### Whittaker's beta ########
print("Whittaker's beta : ")

pdf(paste(workingPath, "figures/all-elevation-brayBin-", otuThr, "-comparison.pdf", sep = ""), width=9, height=5)
attach(mtcars)
par(mfrow=c(1,2), oma=c(4,0,0,0)) 

for (expId in 1:n) {	
    d.brayBin <- whittakerDisList[[ expId ]]
	elevPlotDist <- elevPlotDistList[[ expId ]] 
    
    mantel.brayBin <- mantel(elevPlotDist, d.brayBin, permutations=4999)
   
    # figures
	if (expId == 1) {
		par(mar=c(4,5,2,1))		
		plot(as.vector(elevPlotDist), as.vector(d.brayBin), pch=myshape[expId], col=mypalette[expId], ylim=c(0,1), cex=0.3,
			   xlab="elevation difference (metres)", ylab="Whittaker's beta", main="(a)")	
	} else {
		points(as.vector(elevPlotDist), as.vector(d.brayBin), pch=myshape[expId], col=mypalette[expId], ylim=c(0,1), cex=0.3)
	}
	
	lm.d.brayBin <- lm(as.vector(d.brayBin)~as.vector(elevPlotDist))
	abline(lm.d.brayBin, col=mypalette[expId], cex=3) 
	
	print(paste(matrixNames[expId],"'s mantel statistic r: ", round(mantel.brayBin$statistic, 3), "; significance: ", mantel.brayBin$signif, "; permutations: ", mantel.brayBin$permutations, sep=""))  		
	
}

for (expId in 1:m) {
    d.brayBin <- whittakerNo454DisList[[ expId ]]
	elevPlotDist <- elevPlotDistNo454List[[ expId ]] 

    mantel.brayBin <- mantel(elevPlotDist, d.brayBin, permutations=4999)
    
    # figures
	if (expId == 1) {
		par(mar=c(4,5,2,1))		
		plot(as.vector(elevPlotDist), as.vector(d.brayBin), pch=myshape2[expId], col=mypalette2[expId], ylim=c(0,1), cex=0.3,
			   xlab="elevation difference (metres)", ylab="", main="(b)")	
	} else {
		points(as.vector(elevPlotDist), as.vector(d.brayBin), pch=myshape2[expId], col=mypalette2[expId], ylim=c(0,1), cex=0.3)
	}
	
	lm.d.brayBin <- lm(as.vector(d.brayBin)~as.vector(elevPlotDist))
	abline(lm.d.brayBin, col=mypalette2[expId], cex=3) 
	
	print(paste(matrixNamesNo454[expId],"'s mantel statistic r: ", round(mantel.brayBin$statistic, 3), "; significance: ", mantel.brayBin$signif, "; permutations: ", mantel.brayBin$permutations, sep=""))  		
	
}

par(usr=c(0,1,0,1),xpd=NA)  
legend(-1.38, -0.3, ncol=3, legend=matrixNames, pch=as.numeric(myshape), col=mypalette)
legend(0.1, -0.3, ncol=2, legend=matrixNamesNo454, pch=as.numeric(myshape2), col=mypalette2)
invisible(dev.off())  

####### beta 1 ########
print("beta 1 : ")

pdf(paste(workingPath, "figures/all-elevation-beta1-", otuThr, "-comparison.pdf", sep = ""), width=9, height=5)
attach(mtcars)
par(mfrow=c(1,2), oma=c(4,0,0,0)) 

for (expId in 1:n) {	
    d.beta1 <- beta1DisList[[ expId ]]
	elevPlotDist <- elevPlotDistList[[ expId ]] 
    
    mantel.beta1 <- mantel(elevPlotDist, d.beta1, permutations=4999)
   
    # figures
	if (expId == 1) {
		par(mar=c(4,5,2,1))		
		plot(as.vector(elevPlotDist), as.vector(d.beta1), pch=myshape[expId], col=mypalette[expId], ylim=c(0,1), cex=0.3,
			   xlab="elevation difference (metres)", ylab=expression(paste(""^"1", D[beta], "-1")), main="(a)")	
	} else {
		points(as.vector(elevPlotDist), as.vector(d.beta1), pch=myshape[expId], col=mypalette[expId], ylim=c(0,1), cex=0.3)
	}
	
	lm.d.beta1 <- lm(as.vector(d.beta1)~as.vector(elevPlotDist))
	abline(lm.d.beta1, col=mypalette[expId], cex=3) 
	
	print(paste(matrixNames[expId],"'s mantel statistic r: ", round(mantel.beta1$statistic, 3), "; significance: ", mantel.beta1$signif, "; permutations: ", mantel.beta1$permutations, sep=""))  		
	
}

for (expId in 1:m) {
    d.beta1 <- beta1No454DisList[[ expId ]]
	elevPlotDist <- elevPlotDistNo454List[[ expId ]] 

    #mantel.beta1 <- mantel(elevPlotDist, d.beta1, permutations=4999)
    
    # figures
	if (expId == 1) {
		par(mar=c(4,5,2,1))		
		plot(as.vector(elevPlotDist), as.vector(d.beta1), pch=myshape2[expId], col=mypalette2[expId], ylim=c(0,1), cex=0.3,
			   xlab="elevation difference (metres)", ylab="", main="(b)")	
	} else {
		points(as.vector(elevPlotDist), as.vector(d.beta1), pch=myshape2[expId], col=mypalette2[expId], ylim=c(0,1), cex=0.3)
	}
	
	lm.d.beta1 <- lm(as.vector(d.beta1)~as.vector(elevPlotDist))
	abline(lm.d.beta1, col=mypalette2[expId], cex=3) 
	
	print(paste(matrixNamesNo454[expId],"'s mantel statistic r: ", round(mantel.beta1$statistic, 3), "; significance: ", mantel.beta1$signif, "; permutations: ", mantel.beta1$permutations, sep=""))  		
	
}

par(usr=c(0,1,0,1),xpd=NA)  
legend(-1.38, -0.3, ncol=3, legend=matrixNames, pch=as.numeric(myshape), col=mypalette)
legend(0.1, -0.3, ncol=2, legend=matrixNamesNo454, pch=as.numeric(myshape2), col=mypalette2)
invisible(dev.off())  


####### Horn-Morisita ########
print("Horn-Morisita : ")

pdf(paste(workingPath, "figures/all-elevation-hornMorisita-", otuThr, "-comparison.pdf", sep = ""), width=9, height=5)
attach(mtcars)
par(mfrow=c(1,2), oma=c(4,0,0,0)) 

for (expId in 1:n) {	
    d.hornMorisita <- whittakerDisList[[ expId ]]
	elevPlotDist <- elevPlotDistList[[ expId ]] 
    
    mantel.hornMorisita <- mantel(elevPlotDist, d.hornMorisita, permutations=4999)
   
    # figures
	if (expId == 1) {
		par(mar=c(4,5,2,1))		
		plot(as.vector(elevPlotDist), as.vector(d.hornMorisita), pch=myshape[expId], col=mypalette[expId], ylim=c(0,1), cex=0.3,
			   xlab="elevation difference (metres)", ylab="Horn-Morisita overlap", main="(a)")	
	} else {
		points(as.vector(elevPlotDist), as.vector(d.hornMorisita), pch=myshape[expId], col=mypalette[expId], ylim=c(0,1), cex=0.3)
	}
	
	lm.d.hornMorisita <- lm(as.vector(d.hornMorisita)~as.vector(elevPlotDist))
	abline(lm.d.hornMorisita, col=mypalette[expId], cex=3) 
	
	print(paste(matrixNames[expId],"'s mantel statistic r: ", round(mantel.hornMorisita$statistic, 3), "; significance: ", mantel.hornMorisita$signif, "; permutations: ", mantel.hornMorisita$permutations, sep=""))  		
	
}

for (expId in 1:m) {
    d.hornMorisita <- whittakerNo454DisList[[ expId ]]
	elevPlotDist <- elevPlotDistNo454List[[ expId ]] 

    mantel.hornMorisita <- mantel(elevPlotDist, d.hornMorisita, permutations=4999)
    
    # figures
	if (expId == 1) {
		par(mar=c(4,5,2,1))		
		plot(as.vector(elevPlotDist), as.vector(d.hornMorisita), pch=myshape2[expId], col=mypalette2[expId], ylim=c(0,1), cex=0.3,
			   xlab="elevation difference (metres)", ylab="", main="(b)")	
	} else {
		points(as.vector(elevPlotDist), as.vector(d.hornMorisita), pch=myshape2[expId], col=mypalette2[expId], ylim=c(0,1), cex=0.3)
	}
	
	lm.d.hornMorisita <- lm(as.vector(d.hornMorisita)~as.vector(elevPlotDist))
	abline(lm.d.hornMorisita, col=mypalette2[expId], cex=3) 
	
	print(paste(matrixNamesNo454[expId],"'s mantel statistic r: ", round(mantel.hornMorisita$statistic, 3), "; significance: ", mantel.hornMorisita$signif, "; permutations: ", mantel.hornMorisita$permutations, sep=""))  		
	
}

par(usr=c(0,1,0,1),xpd=NA)  
legend(-1.38, -0.3, ncol=3, legend=matrixNames, pch=as.numeric(myshape), col=mypalette)
legend(0.1, -0.3, ncol=2, legend=matrixNamesNo454, pch=as.numeric(myshape2), col=mypalette2)
invisible(dev.off())  

