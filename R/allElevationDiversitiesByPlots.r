
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
if(!exists("diss.fun")) diss.fun="beta1-1"
if(!exists("taxa.group")) taxa.group="all"


n <- length(matrixNames) 
mypalette <- c("red", "green", "brown", "orange", "blue", "purple", "dark green")
myshape <- seq(0, (0 + n-1))

source("Modules/init.R", local=TRUE)

env.byplot <- getSampleMetaData(TRUE) # by plot

cat("Analysis: elevation difference regression based on", diss.fun, ".\n")

mantel_table <- matrix(0,nrow=n,ncol=4)
rownames(mantel_table) <- c(matrixNames)
colnames(mantel_table) <- c("Mantel statistic $r$", "significance", "R$^2$", "p-value")

####### beta 1 - 1 vs elevation ########
fname <- paste("elev-diff", postfix("all", TRUE, rmSingleton, sep="-"), diss.fun, sep = "-")
pdf(paste(workingPath, figDir, "/", fname, ".pdf", sep = ""), width=9, height=5)
attach(mtcars)
par(mfrow=c(1,2), oma=c(4,0,0,0)) 

for (expId in 1:n) {
  # isPlot, rmSingleton, taxa.group, are fixed in init, when expId == n
  diss <- getDissimilarityMatrix(expId, TRUE, rmSingleton, diss.fun, taxa.group)
  d.diss <- as.dist(diss)
  d.elev.plot <- getElevPlotDist(rownames(diss), env.byplot)
  
  mantel.diss <- mantel(d.elev.plot, d.diss, permutations=4999)
  
  # figures
  if (expId == 1) {
    ######## eDNA #######
    par(mar=c(4,5,2,1))		
    plot(as.vector(d.elev.plot), as.vector(d.diss), pch=myshape[expId], col=mypalette[expId], ylim=c(0,1), cex=0.2,
         xlab="elevation difference (metres)", ylab=expression(paste(""^"1", D[beta], "-1")), main="(a)")	
  } else if (expId == n) {
    ######## vegetation #######
    par(mar=c(4,5,2,1))		
    plot(as.vector(d.elev.plot), as.vector(d.diss), pch=myshape[expId], col=mypalette[expId], ylim=c(0,1), cex=0.2,
         xlab="elevation difference (metres)", ylab="", main="(b)")	
  } else {
    points(as.vector(d.elev.plot), as.vector(d.diss), pch=myshape[expId], col=mypalette[expId], ylim=c(0,1), cex=0.2)
  }
  
  lm.d.diss <- lm(as.vector(d.diss)~as.vector(d.elev.plot))
  abline(lm.d.diss, col=mypalette[expId], cex=3) 
  
  mantel_table[expId,1] <- round(mantel.diss$statistic, 3)
  mantel_table[expId,2] <- mantel.diss$signif
  mantel_table[expId,3] <- formatC(signif(summary(lm.d.diss)$r.squared, digits=3), digits=3,format="fg", flag="#")
  mantel_table[expId,4] <- formatC(signif(summary(lm.d.diss)$coefficients[2,4],digits=3), digits=3,format="fg", flag="#")
  # print(paste(matrixNames[expId],"'s mantel statistic r: ", round(mantel.diss$statistic, 3), "; significance: ", mantel.diss$signif, "; permutations: ", mantel.diss$permutations, sep=""))  		
}

######## combine 2 subfigures #######
par(usr=c(0,1,0,1),xpd=NA)
legend(-1.4, -0.35, ncol=n, legend=matrixNames, pch=as.numeric(myshape), col=mypalette)
invisible(dev.off())  

######## table #######
print(xtable(mantel_table, caption = "Mantel statistic $r$ and their significance using Mantel's test based on 4999 permutations, 
		and R$^2$ and $p$-$value$ for a linear model of the regression of effective $\\beta$ - 1 diversity and difference in elevation.", 
    label = "tab:betaelevationMental", caption.placement = "top"), sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

