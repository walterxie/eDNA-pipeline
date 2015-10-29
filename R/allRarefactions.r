library(RColorBrewer)


if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("isPlot")) isPlot = FALSE # by subplot
if(!exists("otuThr")) otuThr = 97
if(!exists("levels")) levels = rep(c("gamma","alpha","beta"),3)
if(!exists("qs")) qs = rep(0:2,each=3)

subTitles <- c(expression(paste("(e) "^"0", D[gamma])),expression(paste("(a) "^"0", D[alpha])),expression(paste("(c) "^"0", D[beta])),
               expression(paste("(f) "^"1", D[gamma])),expression(paste("(b) "^"1", D[alpha])),expression(paste("(d) "^"1", D[beta])),"(g)","(h)","(i)")

n <- length(matrixNames) 
mypalette <- c("red", "orange", "green", "purple", "blue", "brown")
myshape <- seq(15, (15 + n-2))

######## rarefaction #######
pdf(paste(workingPath, figDir, "/", postfix("edna-rarefaction", isPlot, rmSingleton, sep="-"), ".pdf", sep = ""), width=6, height=9)
attach(mtcars)
par(mfrow=c(3,2), oma=c(3,0,0,0))

for (dFId in c(2,5,3,6,1,4)) { #length(levels)
  maxX <- 0
  maxY <- 0
  
  for (expId in 1:(n-1)) {
    rarefactionTable <- getRarefactionTable(expId, isPlot, rmSingleton)
    # remove prefix size.
    sampleSizesSeq <- gsub("^.*?\\.","",colnames(rarefactionTable)) 
    
    maxX <- max(maxX, max(as.numeric(sampleSizesSeq)))
    maxY <- max(maxY, max(as.numeric(rarefactionTable[dFId,])))
  }
  
  for (expId in 1:(n-1)) {
    rarefactionTable <- getRarefactionTable(expId, isPlot, rmSingleton)
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
      
      #pdf(paste("figures/", postfix("all-rarefaction", isPlot, rmSingleton, sep="-"), "-", paste(levels, qs, sep="")[dFId], ".pdf", sep = ""))
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
legend(-1.1, -0.3, ncol=6, legend=matrixNames[1:(n-1)], col=mypalette, pch=as.numeric(myshape))

invisible(dev.off())