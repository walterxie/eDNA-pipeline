library(RColorBrewer)

# change config below
#figDir <- "figures"
#sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
#setwd(sourcePath)
#workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
#matrixNames <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name   
#levels = rep(c("gamma","alpha","beta"),3)
#qs = rep(0:2,each=3)

if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")
if(!exists("levels")) stop("levels of Jost diversity are missing !")
if(!exists("qs")) stop("qs of Jost diversity are missing !")
if(!exists("rmSingleton")) stop("rmSingleton flag is missing !")

subTitles <- c(expression(paste("(e) "^"0", D[gamma])),expression(paste("(a) "^"0", D[alpha])),expression(paste("(c) "^"0", D[beta])),
               expression(paste("(f) "^"1", D[gamma])),expression(paste("(b) "^"1", D[alpha])),expression(paste("(d) "^"1", D[beta])),"(g)","(h)","(i)")

n <- length(matrixNames) 
#mypalette <- rainbow(n)
mypalette <- c("red", "orange", "green", "purple", "blue", "brown")
myshape <- seq(0, (0 + n-1))
	
stringBySubOrPlot <- "-by-subplot"
otuThrSeq <- 90:100

pdf(paste(workingPath, figDir, "/all-diver-otu-log.pdf", sep = ""), width=6, height=9)
attach(mtcars)
par(mfrow=c(3,2), oma=c(3,0,0,0))

for (dFId in c(2,5,3,6,1,4)) { #length(levels)
    rarefractionTable <- NULL    
    maxY <- 0 
    minY <- 999999 
    for (expId in 1:n) {
		matrixName <- paste(matrixNames[expId], stringBySubOrPlot, sep = "")

		if (rmSingleton) 
			matrixName <- paste("nonsingleton", matrixName, sep = "-")

		inputRDT <- paste(workingPath, matrixNames[expId], "/", matrixName, "-otus-thre-table.csv", sep = "")
		rarefractionTable <- read.csv(file=inputRDT, head=TRUE, sep=",", row.names=paste(levels, qs, sep=""), check.names=FALSE)
		
		maxY <- max(maxY, max(as.numeric(rarefractionTable[dFId,])))
		minY <- min(minY, min(as.numeric(rarefractionTable[dFId,])))
	}
	
	rarefractionTable <- NULL 	
	for (expId in 1:n) {
		matrixName <- paste(matrixNames[expId], stringBySubOrPlot, sep = "")
		
		if (rmSingleton) 
			matrixName <- paste("nonsingleton", matrixName, sep = "-")

    	inputRDT <- paste(workingPath, matrixNames[expId], "/", matrixName, "-otus-thre-table.csv", sep = "")
		rarefractionTable <- read.csv(file=inputRDT, head=TRUE, sep=",", row.names=paste(levels, qs, sep=""), check.names=FALSE)        
	
		if (expId == 1) {
			if (levels[dFId] == "gamma") {
			xlab="OTUs % threshold"
			par(mar=c(4,4,3,2)) 
			} else {
			xlab=""
			par(mar=c(4,4,3,2))
			}
		
			#pdf(paste("figures/diver-otu-", paste(levels, qs, sep="")[dFId], ".pdf", sep = ""))
			if (levels[dFId] == "beta") {
				plot(otuThrSeq, rarefractionTable[dFId,], pch=myshape[expId], col=mypalette[expId], xlim=rev(range(otuThrSeq)), ylim=c(0,maxY),
					   xlab=xlab, ylab=paste(levels, "-", qs, " diversity", sep="")[dFId], main=subTitles[dFId]) 				  
			} else {
				plot(otuThrSeq, log10(rarefractionTable[dFId,]), pch=myshape[expId], col=mypalette[expId], xlim=rev(range(otuThrSeq)), ylim=c(log10(minY),log10(maxY)), yaxt="n", 
					   xlab=xlab, ylab=paste(levels, "-", qs, " diversity", sep="")[dFId], main=subTitles[dFId]) 	
					   
				aty <- axTicks(2)
				labels <- sapply(aty,function(i) round(10^i, 0))
				axis(2,at=aty,labels=labels)			  			
			}
			
		} else {
			if (levels[dFId] == "beta") {
				points(otuThrSeq, rarefractionTable[dFId,], pch=myshape[expId], col=mypalette[expId]) 	  
			} else {
				points(otuThrSeq, log10(rarefractionTable[dFId,]), pch=myshape[expId], col=mypalette[expId]) 			  			
			}
			
		}
	
		lty=2
		if (expId < 3 | expId == 5) { # 16S, 18S and COI solid 
		   lty=1
		}
	
		if (levels[dFId] == "beta") {
			lines(otuThrSeq, rarefractionTable[dFId,], lty=lty, col=mypalette[expId]) 
		} else {
			lines(otuThrSeq, log10(rarefractionTable[dFId,]), lty=lty, col=mypalette[expId]) 			  			
		}			
		
		abline(v=97, lty=2, col="gray60")
	}    
} 

par(usr=c(0,1,0,1),xpd=NA)
legend(-1.2, -0.3, ncol=6, legend=matrixNames, col=mypalette, pch=as.numeric(myshape))

invisible(dev.off())
