library(reshape2)

if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97

n <- length(matrixNames) 

source("Modules/init.R", local=TRUE)

######## within plots vs between plots #######	
for (expId in 1:(n-1)) {	
  beta1_1 <- getBeta1Minus1(expId, FALSE, rmSingleton)
    
  beta1_1.melt <- melt(as.matrix(beta1_1))
  beta1_1.melt <- beta1_1.melt[beta1_1.melt$value>0,]
  
  beta1_1.melt$dataset <- matrixNames[expId]
  beta1_1.melt$comp <- "between" 
  # plot-subplot
  beta1_1.melt[,1] %in% beta1_1.melt[,2]
  beta1_1.melt[, "comp"] <- "within"

  
  if (expId==1) {
    diver <- beta1_1.melt
  } else {
    diver <- rbind(diver, beta1_1.melt)
  }
  
}

####### within plots vs between plots #######

diverInAPlot$comp <- paste(diverInAPlot$region, "-within", sep = "") 
diverBetweenPlot$comp <- paste(diverBetweenPlot$region, "-between", sep = "")
diver <- rbind(diverInAPlot, diverBetweenPlot)
comparasion <- paste(rep(matrixNames, each = 2), rep(c("-within", "-between"), times = length(matrixNames)), sep = "")
diver$comp = factor(diver$comp,comparasion)

pdf(paste(workingPath, figDir, "/all-beta1-", otuThr, "-comparison.pdf", sep = ""), width=9, height=6)

par(mar=c(3,5,1,1)) 
boxplot(beta1-1 ~ comp, data=diver, at = sort(c(seq(1,(n-1)*3+1,by=3),seq(2,(n-1)*3+2,by=3))), xaxt = 'n',
    col=rep(c("red","blue"),times=n), ylab =expression(paste(""^"1", D[beta], "-1")))
axis(1, at=seq(1.5,(n-1)*3+1.5,by=3), labels=matrixNames)
legend("bottomright", legend=c("Within","Between"), fill=c("red","blue"))
		
invisible(dev.off())



