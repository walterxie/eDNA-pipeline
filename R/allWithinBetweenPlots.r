library(reshape2)

if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("diss.fun")) diss.fun="beta1-1"
if(!exists("taxa.group")) taxa.group="all"

n <- length(matrixNames) 

source("Modules/init.R", local=TRUE)

cat("Analysis:", taxa.group, "taxa group", diss.fun, "dissimilarity within plots vs between plots. \n")

######## within plots vs between plots #######	
for (expId in 1:(n-1)) {	
  diss <- getDissimilarityMatrix(expId, FALSE, rmSingleton, diss.fun, taxa.group)
    
  diss.melt <- melt(as.matrix(diss))
  diss.melt <- diss.melt[diss.melt$value>0,]
  
  diss.melt$dataset <- matrixNames[expId]
  diss.melt$comp <- "between" 
  # plot-subplot
  diss.melt[getPlot(diss.melt[,1])==getPlot(diss.melt[,2]), "comp"] <- "within"

  if (expId==1) {
    diver <- diss.melt
  } else {
    diver <- rbind(diver, diss.melt)
  }
}

####### within plots vs between plots #######
diver$comp <- paste(diver$dataset, diver$comp, sep = "-")
diver$comp = factor(diver$comp, unique(diver$comp))
n_edna=n-1

fname <- paste("within-between", postfix(taxa.group, FALSE, rmSingleton, sep="-"), diss.fun, sep = "-")
pdf(paste(workingPath, figDir, "/", fname, ".pdf", sep = ""), width=9, height=6)

par(mar=c(3,5,1,1)) 
boxplot(value ~ comp, data=diver, at = sort(c(seq(1,(n_edna-1)*3+1,by=3),seq(2,(n_edna-1)*3+2,by=3))), xaxt = 'n',
    col=rep(c("red","blue"),times=n_edna), ylab =expression(paste(""^"1", D[beta], "-1")))
axis(1, at=seq(1.5,(n_edna-1)*3+1.5,by=3), labels=matrixNames[1:n_edna])
legend("topleft", legend=c("Within","Between"), fill=c("red","blue"))
		
invisible(dev.off())



