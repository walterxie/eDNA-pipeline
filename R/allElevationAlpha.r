# "-by-plot" trigger merge 2 subplots columns
library(RColorBrewer)
library(vegan)
library(vegetarian)
library(xtable)

alpha1 <- function(communityMatrix) {    
	####### alpha 1 ########
	# including diagonal
    d.alpha1 <- rep(0,nrow(communityMatrix))	
	for(i in 1:nrow(communityMatrix)){				
		d.alpha1[i] <- d(communityMatrix[i,],lev="gamma",q=1)				
	}
	
	return (d.alpha1)
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
alpha1DisList <- list()
elevPlotDistList <- list()

max454 = 0
for (expId in 1:n) {	
    communityMatrix <- init(expId, otuThr, "-by-subplot")
    elevPlotDist <- getElevSample(rownames(communityMatrix))
    
    elevPlotDistList[[ expId ]] <- elevPlotDist	
    
	####### beta 1 ########
	# (45 pairs for 10 plots)
    d.alpha1 <- alpha1(communityMatrix)
    #mantel.alpha1 <- mantel(elevPlotDist, d.alpha1, permutations=4999)
    
    alpha1DisList[[ expId ]] <- d.alpha1

	max.tmp <- max(d.alpha1)
 	max454 <- max(max454, max.tmp)

    print(paste(matrixNames[expId], " (", otuThr, "%) : total = ", sum(communityMatrix), "; min sample per subplot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))
}

######## non 454 #######	
alpha1No454DisList <- list()
elevPlotDistNo454List <- list()

maxnon454 = 0
for (expId in 1:m) {    	
    communityMatrix <- initNon454ByPlot(expId, otuThr)
    elevPlotDist <- getElevSample(rownames(communityMatrix))
    
    elevPlotDistNo454List[[ expId ]] <- elevPlotDist	
    	
	####### alpha 1 ########
	# (45 pairs for 10 plots)
    d.alpha1 <- alpha1(communityMatrix)
    #mantel.alpha1 <- mantel(elevPlotDist, d.alpha1, permutations=4999)
    
    alpha1No454DisList[[ expId ]] <- d.alpha1
    
 	max.tmp <- max(d.alpha1)
 	maxnon454 <- max(maxnon454, max.tmp)
    
    print(paste(matrixNamesNo454[expId], " (", otuThr, "%) : total = ", sum(communityMatrix), "; min sample per subplot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))
}

####### alpha 1 ########
print("alpha 1 : ")
sig_table <- matrix(0,nrow=(n+m),ncol=2)
rownames(sig_table) <- c(matrixNames, matrixNamesNo454)
colnames(sig_table) <- c("R$^2$", "p-value")

pdf(paste(workingPath, "figures/all-elevation-alpha1-", otuThr, "-comparison.pdf", sep = ""), width=9, height=5)
attach(mtcars)
par(mfrow=c(1,2), oma=c(4,0,0,0)) 

for (expId in 1:n) {	
	d.alpha1 <- alpha1DisList[[ expId ]]
	elevPlotDist <- elevPlotDistList[[ expId ]] 

	# figures
	if (expId == 1) {
		par(mar=c(4,5,2,1))		
		plot(elevPlotDist[,1], log10(as.vector(d.alpha1)), pch=myshape[expId], col=mypalette[expId], ylim=c(0,log10(3000)), cex=0.8,
			   xlab="elevation (metres)", ylab=expression(paste(""^"1", D[alpha])), main="(a)", yaxt="n")	

		aty <- axTicks(2)
		ylabels <- sapply(aty,function(i) 10^i)
		axis(2,at=aty,labels=ylabels)    

	} else {
		points(elevPlotDist[,1], log10(as.vector(d.alpha1)), pch=myshape[expId], col=mypalette[expId], ylim=c(0,log10(3000)), cex=0.8)
	}

	lm.d.alpha1 <- lm(log10(as.vector(d.alpha1))~elevPlotDist[,1])
	abline(lm.d.alpha1, col=mypalette[expId], cex=3) 

	sig_table[expId,1] <- formatC(signif(summary(lm.d.alpha1)$r.squared, digits=3), digits=3,format="fg", flag="#")
	sig_table[expId,2] <- formatC(signif(summary(lm.d.alpha1)$coefficients[2,4],digits=3), digits=3,format="fg", flag="#")
}

for (expId in 1:m) {
	d.alpha1 <- alpha1No454DisList[[ expId ]]
	elevPlotDist <- elevPlotDistNo454List[[ expId ]] 

	# figures
	if (expId == 1) {
		par(mar=c(4,5,2,1))		
		plot(elevPlotDist[,1], log10(as.vector(d.alpha1)), pch=myshape2[expId], col=mypalette2[expId], ylim=c(0,log10(100)), cex=0.8,
			   xlab="elevation (metres)", ylab="", main="(b)", yaxt="n")	
						
		aty <- axTicks(2)
		ylabels <- sapply(aty,function(i) 10^i)
		axis(2,at=aty,labels=ylabels)    

	} else {
		points(elevPlotDist[,1], log10(as.vector(d.alpha1)), pch=myshape2[expId], col=mypalette2[expId], ylim=c(0,log10(100)), cex=0.8)
	}

	lm.d.alpha1 <- lm(log10(as.vector(d.alpha1))~elevPlotDist[,1])
	abline(lm.d.alpha1, col=mypalette2[expId], cex=3) 

	sig_table[expId+n,1] <- formatC(signif(summary(lm.d.alpha1)$r.squared, digits=3), digits=3,format="fg", flag="#")
	sig_table[expId+n,2] <- formatC(signif(summary(lm.d.alpha1)$coefficients[2,4],digits=3), digits=3,format="fg", flag="#")

}

par(usr=c(0,1,0,1),xpd=NA)  
legend(-1.38, -0.3, ncol=3, legend=matrixNames, pch=as.numeric(myshape), col=mypalette)
legend(0.1, -0.3, ncol=2, legend=matrixNamesNo454, pch=as.numeric(myshape2), col=mypalette2)
invisible(dev.off())  

print(xtable(sig_table),sanitize.text.function=function(x){x})


