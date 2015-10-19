
library(RColorBrewer)
library(vegan)
library(vegetarian)
library(xtable)


if(!exists("tableFile")) stop("table file is missing !")
if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97

n <- length(matrixNames) 
mypalette <- c("red", "orange", "green", "purple", "blue", "brown", "black")
myshape <- seq(0, (0 + n-1))

source("Modules/init.R", local=TRUE)

env <- getSampleMetaData(TRUE) # by plot

####### beta 1 ########
print("beta 1 : ")
mantel_table <- matrix(0,nrow=(n+m),ncol=4)
rownames(mantel_table) <- c(matrixNames, matrixNamesNo454)
colnames(mantel_table) <- c("Mantel statistic $r$", "significance", "R$^2$", "p-value")

pdf(paste(workingPath, figDir, "/all-elev-diff-beta1-1.pdf", sep = ""), width=9, height=5)
attach(mtcars)
par(mfrow=c(1,2), oma=c(4,0,0,0)) 

######## eDNA #######
for (expId in 1:(n-1)) {	
    d.beta1_1 <- as.dist(getBeta1Minus1(expId, TRUE, rmSingleton))
	d.elev.plot <- elevPlotDistList[[ expId ]] 
    
    mantel.beta1_1 <- mantel(elevPlotDist, d.beta1_1, permutations=4999)
   
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

######## vegetation #######
d.beta1_1 <- as.dist(getBeta1Minus1(n, TRUE, FALSE))
d.elev.plot <- elevPlotDistList[[ expId ]] 

mantel.beta1_1 <- mantel(elevPlotDist, d.beta1_1, permutations=4999)
    
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
	


par(usr=c(0,1,0,1),xpd=NA)  
legend(-1.38, -0.3, ncol=3, legend=matrixNames, pch=as.numeric(myshape), col=mypalette)
legend(0.1, -0.3, ncol=2, legend=matrixNamesNo454, pch=as.numeric(myshape2), col=mypalette2)
invisible(dev.off())  

print(xtable(mantel_table, caption = "Mantel statistic $r$ and their significance using Mantel's test based on 4999 permutations, 
		and R$^2$ and $p$-$value$ for a linear model of the regression of effective $\\beta$ diversity and difference in elevation.", 
		label = "tab:betaelevationMental", caption.placement = "top"), sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

