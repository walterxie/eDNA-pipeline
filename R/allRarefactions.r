# change config below
workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
experiments <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name
matrixNames <-  experiments
subTitles <- c(expression(paste("(e) "^"0", D[gamma])),expression(paste("(a) "^"0", D[alpha])),expression(paste("(c) "^"0", D[beta])),
               expression(paste("(f) "^"1", D[gamma])),expression(paste("(b) "^"1", D[alpha])),expression(paste("(d) "^"1", D[beta])),"(g)","(h)","(i)")

n <- length(matrixNames) 
mypalette <- c("red", "orange", "green", "purple", "blue", "brown")
myshape <- seq(15, (15 + n -1))

levels = rep(c("gamma","alpha","beta"),3)
qs = rep(0:2,each=3)
otuThr = 97

stringBySubOrPlot <- "-by-subplot"

pdf(paste(workingPath, "figures/all-rarefaction-", otuThr, ".pdf", sep = ""), width=6, height=9)
attach(mtcars)
par(mfrow=c(3,2), oma=c(3,0,0,0))

for (dFId in c(2,5,3,6,1,4)) { #length(levels)
    maxX <- 0
    maxY <- 0
    
    for (expId in 1:length(matrixNames)) {
        matrixName <- paste(matrixNames[expId], stringBySubOrPlot, sep = "")
        inputRDT <- paste(workingPath, experiments[expId], "/", matrixName, "-", otuThr, "-rarefaction-table.csv", sep = "")
	
		rarefactionTable <- read.csv(file=inputRDT, head=TRUE, sep=",", row.names=paste(levels, qs, sep=""), check.names=FALSE)
		 # remove prefix size.
	    sampleSizesSeq <- gsub("^.*?\\.","",colnames(rarefactionTable)) 
	    
		maxX <- max(maxX, max(as.numeric(sampleSizesSeq)))
        maxY <- max(maxY, max(as.numeric(rarefactionTable[dFId,])))
    }
        
	for (expId in 1:length(matrixNames)) {
        matrixName <- paste(matrixNames[expId], stringBySubOrPlot, sep = "")
        inputRDT <- paste(workingPath, experiments[expId], "/", matrixName, "-", otuThr, "-rarefaction-table.csv", sep = "")
	
		rarefactionTable <- read.csv(file=inputRDT, head=TRUE, sep=",", row.names=paste(levels, qs, sep=""), check.names=FALSE)
	    # remove prefix size.
	    sampleSizesSeq <- gsub("^.*?\\.","",colnames(rarefactionTable)) 
	
		if (expId == 1) {
			if (levels[dFId] == "gamma") {
				xlab="Sample size per site"
				par(mar=c(4,4,3,2)) 
			} else {
				xlab=""
				par(mar=c(4,4,3,2))
			}
			
			#pdf(paste("figures/all-rarefaction-", otuThr, "-", paste(levels, qs, sep="")[dFId], ".pdf", sep = ""))
	        par(mar=c(4,4,3,2))
			plot(sampleSizesSeq, rarefactionTable[dFId,], pch=myshape[expId], col=mypalette[expId], xlim=c(0,maxX), ylim=c(0,maxY), 
				   xlab=xlab, ylab=paste(levels, "-", qs, " diversity", sep="")[dFId], main=subTitles[dFId]) 				  
			
		} else {
			points(sampleSizesSeq, rarefactionTable[dFId,], pch=myshape[expId], col=mypalette[expId]) 
		}
	
		lines(sampleSizesSeq, rarefactionTable[dFId,], lty=2, col=mypalette[expId]) 
	}    

} 

par(usr=c(0,1,0,1),xpd=NA)
legend(-1.2, -0.3, ncol=6, legend=matrixNames, col=mypalette, pch=as.numeric(myshape))

invisible(dev.off())