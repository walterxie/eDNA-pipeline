library(reshape2)

if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97

n <- length(matrixNames) 

source("Modules/init.R", local=TRUE)

cat("Analysis: beta1 - 1 of within plots vs between plots. \n")

######## within plots vs between plots #######	
for (expId in 1:(n-1)) {	
  beta1_1 <- getBeta1Minus1(expId, FALSE, rmSingleton)
    
  beta1_1.melt <- melt(as.matrix(beta1_1))
  beta1_1.melt <- beta1_1.melt[beta1_1.melt$value>0,]
  
  beta1_1.melt$dataset <- matrixNames[expId]
  beta1_1.melt$comp <- "between" 
  # plot-subplot
  beta1_1.melt[getPlot(beta1_1.melt[,1])==getPlot(beta1_1.melt[,2]), "comp"] <- "within"

  if (expId==1) {
    diver <- beta1_1.melt
  } else {
    diver <- rbind(diver, beta1_1.melt)
  }
}

####### within plots vs between plots #######
diver$dataset = factor(diver$dataset, unique(diver$dataset))
diver$comp = factor(diver$comp, unique(diver$comp))
n_edna=n-1

pdf(paste(workingPath, figDir, "/", postfix("all-within-between", FALSE, rmSingleton, sep="-"), "-beta1-1.pdf", sep = ""), width=9, height=6)

par(mar=c(3,5,1,1)) 
boxplot(value ~ interaction(dataset, comp), data=diver, at = sort(c(seq(1,(n_edna-1)*3+1,by=3),seq(2,(n_edna-1)*3+2,by=3))), xaxt = 'n',
    col=rep(c("red","blue"),times=n_edna), ylab =expression(paste(""^"1", D[beta], "-1")))
axis(1, at=seq(1.5,(n_edna-1)*3+1.5,by=3), labels=matrixNames[1:n_edna])
legend("topleft", legend=c("Within","Between"), fill=c("red","blue"))
		
invisible(dev.off())



