# "-by-plot" trigger merge 2 subplots columns
library(RColorBrewer)
library(vegan)
library(vegetarian)
library(xtable)

# change config below
#figDir <- "figures"
#sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
#setwd(sourcePath)
#workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
#matrixNames <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name   
#matrixNamesNo454 <-  c("seedlings","trees","invertebrates","birds") # need expId correct for "birds","seedlings" 
#levels = rep(c("gamma","alpha","beta"),3)
#qs = rep(0:2,each=3)

if(!exists("tableFile")) stop("table file is missing !")
if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")
if(!exists("matrixNamesNo454")) stop("matrix names of traditional methods are missing !")
if(!exists("levels")) stop("levels of Jost diversity are missing !")
if(!exists("qs")) stop("qs of Jost diversity are missing !")
if(!exists("rmSingleton")) stop("rmSingleton flag is missing !")

n <- length(matrixNames) 
mypalette <- c("red", "orange", "green", "purple", "blue", "brown")
myshape <- seq(0, (0 + n-1))

m <- length(matrixNamesNo454)
mypalette2 <- brewer.pal(m,"Dark2")
myshape2 <- seq(0, (0 + m-1))

otuThr = 97

source("Modules/init.R", local=TRUE)

######## 454 #######
jaccardDisList <- list()
beta1DisList <- list()
hornMorisitaDisList <- list()
elevPlotDistList <- list()

for (expId in 1:n) {	
    communityMatrix <- init(expId, otuThr, "-by-plot")
    elevPlotDist <- getElevSampleDist(rownames(communityMatrix)) 
    
	######## jaccard #######	
	d.jaccard <- vegdist(communityMatrix, method="jaccard")
	#mantel.jaccard <- mantel(elevPlotDist, d.jaccard, permutations=4999)
    
    jaccardDisList[[ expId ]] <- d.jaccard
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
jaccardNo454DisList <- list()
beta1No454DisList <- list()
hornMorisitaNo454DisList <- list()
elevPlotDistNo454List <- list()

for (expId in 1:m) {    	
    communityMatrix <- initNon454ByPlot(expId, otuThr)
    elevPlotDist <- getElevSampleDist(rownames(communityMatrix)) 
    
	######## jaccard #######	
	d.jaccard <- vegdist(communityMatrix, method="jaccard")
	#mantel.jaccard <- mantel(elevPlotDist, d.jaccard, permutations=4999)
    
    jaccardNo454DisList[[ expId ]] <- d.jaccard
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


####### beta 1 ########
print("beta 1 : ")
mantel_table <- matrix(0,nrow=(n+m),ncol=4)
rownames(mantel_table) <- c(matrixNames, matrixNamesNo454)
colnames(mantel_table) <- c("Mantel statistic $r$", "significance", "R$^2$", "p-value")

pdf(paste(workingPath, figDir, "/all-elevation-beta1-", otuThr, "-comparison.pdf", sep = ""), width=9, height=5)
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
		
	mantel_table[expId,1] <- round(mantel.beta1$statistic, 3)
	mantel_table[expId,2] <- mantel.beta1$signif
	mantel_table[expId,3] <- formatC(signif(summary(lm.d.beta1)$r.squared, digits=3), digits=3,format="fg", flag="#")
	mantel_table[expId,4] <- formatC(signif(summary(lm.d.beta1)$coefficients[2,4],digits=3), digits=3,format="fg", flag="#")
#    print(paste(matrixNames[expId],"'s mantel statistic r: ", round(mantel.beta1$statistic, 3), "; significance: ", mantel.beta1$signif, "; permutations: ", mantel.beta1$permutations, sep=""))  		
	
}

for (expId in 1:m) {
    d.beta1 <- beta1No454DisList[[ expId ]]
	elevPlotDist <- elevPlotDistNo454List[[ expId ]] 

    mantel.beta1 <- mantel(elevPlotDist, d.beta1, permutations=4999)
    
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
	
	mantel_table[(expId+n),1] <- round(mantel.beta1$statistic, 3)
	mantel_table[(expId+n),2] <- mantel.beta1$signif
	mantel_table[(expId+n),3] <- formatC(signif(summary(lm.d.beta1)$r.squared, digits=3), digits=3,format="fg", flag="#")
	mantel_table[(expId+n),4] <- formatC(signif(summary(lm.d.beta1)$coefficients[2,4],digits=3), digits=3,format="fg", flag="#")
#	print(paste(matrixNamesNo454[expId],"'s mantel statistic r: ", round(mantel.beta1$statistic, 3), "; significance: ", mantel.beta1$signif, "; permutations: ", mantel.beta1$permutations, sep=""))  		
	
}

par(usr=c(0,1,0,1),xpd=NA)  
legend(-1.38, -0.3, ncol=3, legend=matrixNames, pch=as.numeric(myshape), col=mypalette)
legend(0.1, -0.3, ncol=2, legend=matrixNamesNo454, pch=as.numeric(myshape2), col=mypalette2)
invisible(dev.off())  

print(xtable(mantel_table, caption = "Mantel statistic $r$ and their significance using Mantel's test based on 4999 permutations, 
		and R$^2$ and $p$-$value$ for a linear model of the regression of effective $\\beta$ diversity and difference in elevation.", 
		label = "tab:betaelevationMental", caption.placement = "top"), sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

####### jaccard ########
print("jaccard : ")
mantel_table <- matrix(0,nrow=(n+m),ncol=4)
rownames(mantel_table) <- c(matrixNames, matrixNamesNo454)
colnames(mantel_table) <- c("Mantel statistic $r$", "significance", "R$^2$", "p-value")

pdf(paste(workingPath, figDir, "/all-elevation-jaccard-", otuThr, "-comparison.pdf", sep = ""), width=9, height=5)
attach(mtcars)
par(mfrow=c(1,2), oma=c(4,0,0,0)) 

for (expId in 1:n) {	
    d.jaccard <- jaccardDisList[[ expId ]]
	elevPlotDist <- elevPlotDistList[[ expId ]] 
    
    mantel.jaccard <- mantel(elevPlotDist, d.jaccard, permutations=4999)
   
    # figures
	if (expId == 1) {
		par(mar=c(4,5,2,1))		
		plot(as.vector(elevPlotDist), as.vector(d.jaccard), pch=myshape[expId], col=mypalette[expId], ylim=c(0,1), cex=0.3,
			   xlab="elevation difference (metres)", ylab="Jaccard index", main="(a)")	
	} else {
		points(as.vector(elevPlotDist), as.vector(d.jaccard), pch=myshape[expId], col=mypalette[expId], ylim=c(0,1), cex=0.3)
	}
	
	lm.d.jaccard <- lm(as.vector(d.jaccard)~as.vector(elevPlotDist))
	abline(lm.d.jaccard, col=mypalette[expId], cex=3) 
	
	mantel_table[expId,1] <- round(mantel.jaccard$statistic, 3)
	mantel_table[expId,2] <- mantel.jaccard$signif
	mantel_table[expId,3] <- formatC(signif(summary(lm.d.jaccard)$r.squared, digits=3), digits=3,format="fg", flag="#")
	mantel_table[expId,4] <- formatC(signif(summary(lm.d.jaccard)$coefficients[2,4],digits=3), digits=3,format="fg", flag="#")
#	print(paste(matrixNames[expId],"'s mantel statistic r: ", round(mantel.jaccard$statistic, 3), "; significance: ", mantel.jaccard$signif, "; permutations: ", mantel.jaccard$permutations, sep=""))  		
	
}

for (expId in 1:m) {
    d.jaccard <- jaccardNo454DisList[[ expId ]]
	elevPlotDist <- elevPlotDistNo454List[[ expId ]] 

    mantel.jaccard <- mantel(elevPlotDist, d.jaccard, permutations=4999)
    
    # figures
	if (expId == 1) {
		par(mar=c(4,5,2,1))		
		plot(as.vector(elevPlotDist), as.vector(d.jaccard), pch=myshape2[expId], col=mypalette2[expId], ylim=c(0,1), cex=0.3,
			   xlab="elevation difference (metres)", ylab="", main="(b)")	
	} else {
		points(as.vector(elevPlotDist), as.vector(d.jaccard), pch=myshape2[expId], col=mypalette2[expId], ylim=c(0,1), cex=0.3)
	}
	
	lm.d.jaccard <- lm(as.vector(d.jaccard)~as.vector(elevPlotDist))
	abline(lm.d.jaccard, col=mypalette2[expId], cex=3) 
	
	mantel_table[(expId+n),1] <- round(mantel.jaccard$statistic, 3)
	mantel_table[(expId+n),2] <- mantel.jaccard$signif
	mantel_table[(expId+n),3] <- formatC(signif(summary(lm.d.jaccard)$r.squared, digits=3), digits=3,format="fg", flag="#")
	mantel_table[(expId+n),4] <- formatC(signif(summary(lm.d.jaccard)$coefficients[2,4],digits=3), digits=3,format="fg", flag="#")
#	print(paste(matrixNamesNo454[expId],"'s mantel statistic r: ", round(mantel.jaccard$statistic, 3), "; significance: ", mantel.jaccard$signif, "; permutations: ", mantel.jaccard$permutations, sep=""))  		
	
}

par(usr=c(0,1,0,1),xpd=NA)  
legend(-1.38, -0.3, ncol=3, legend=matrixNames, pch=as.numeric(myshape), col=mypalette)
legend(0.1, -0.3, ncol=2, legend=matrixNamesNo454, pch=as.numeric(myshape2), col=mypalette2)
invisible(dev.off())  

print(xtable(mantel_table, caption = "Mantel statistic $r$ and their significance using Mantel's test based on 4999 permutations, 
		and R$^2$ and $p$-$value$ for a linear model of the regression of difference in community measured by Jaccard and difference in elevation.", 
		label = "tab:JaccardElevationMental", caption.placement = "top"), sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

####### Horn-Morisita ########
print("Horn-Morisita : ")
mantel_table <- matrix(0,nrow=(n+m),ncol=4)
rownames(mantel_table) <- c(matrixNames, matrixNamesNo454)
colnames(mantel_table) <- c("Mantel statistic $r$", "significance", "R$^2$", "p-value")

pdf(paste(workingPath, figDir, "/all-elevation-hornMorisita-", otuThr, "-comparison.pdf", sep = ""), width=9, height=5)
attach(mtcars)
par(mfrow=c(1,2), oma=c(4,0,0,0)) 

for (expId in 1:n) {	
    d.hornMorisita <- hornMorisitaDisList[[ expId ]]
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
	
	mantel_table[expId,1] <- round(mantel.hornMorisita$statistic, 3)
	mantel_table[expId,2] <- mantel.hornMorisita$signif
	mantel_table[expId,3] <- formatC(signif(summary(lm.d.hornMorisita)$r.squared, digits=3), digits=3,format="fg", flag="#")
	mantel_table[expId,4] <- formatC(signif(summary(lm.d.hornMorisita)$coefficients[2,4],digits=3), digits=3,format="fg", flag="#")
#	print(paste(matrixNames[expId],"'s mantel statistic r: ", round(mantel.hornMorisita$statistic, 3), "; significance: ", mantel.hornMorisita$signif, "; permutations: ", mantel.hornMorisita$permutations, sep=""))  		
	
}

for (expId in 1:m) {
    d.hornMorisita <- hornMorisitaNo454DisList[[ expId ]]
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
	
	mantel_table[(expId+n),1] <- round(mantel.hornMorisita$statistic, 3)
	mantel_table[(expId+n),2] <- mantel.hornMorisita$signif
	mantel_table[(expId+n),3] <- formatC(signif(summary(lm.d.hornMorisita)$r.squared, digits=3), digits=3,format="fg", flag="#")
	mantel_table[(expId+n),4] <- formatC(signif(summary(lm.d.hornMorisita)$coefficients[2,4],digits=3), digits=3,format="fg", flag="#")
#	print(paste(matrixNamesNo454[expId],"'s mantel statistic r: ", round(mantel.hornMorisita$statistic, 3), "; significance: ", mantel.hornMorisita$signif, "; permutations: ", mantel.hornMorisita$permutations, sep=""))  		
	
}

par(usr=c(0,1,0,1),xpd=NA)  
legend(-1.38, -0.3, ncol=3, legend=matrixNames, pch=as.numeric(myshape), col=mypalette)
legend(0.1, -0.3, ncol=2, legend=matrixNamesNo454, pch=as.numeric(myshape2), col=mypalette2)
invisible(dev.off())  

print(xtable(mantel_table, caption = "Mantel statistic $r$ and their significance using Mantel's test based on 4999 permutations, 
		and R$^2$ and $p$-$value$ for a linear model of the regression of difference in community measured by Horn-Morisita and difference in elevation.", 
		label = "tab:HornMorisitaElevationMental", caption.placement = "top"), sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

