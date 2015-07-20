library(vegetarian)

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

n <- length(matrixNames) 

otuThr = 97

source("Modules/init.R", local=TRUE)

######## within plots vs between plots #######	
for (expId in 1:n) {	
    # "-by-plot" trigger merge 2 subplots columns
    communityMatrix <- init(expId, otuThr, "-by-subplot")
    
    # by subplots
    rowsNum <- nrow(communityMatrix)		
	d.subplot<-matrix(0,nrow=rowsNum,ncol=length(levels))
	for(i in 1:rowsNum){
		# levels = [ "gamma" "alpha" "beta"  "gamma" "alpha" "beta"  "gamma" "alpha" "beta" ]
		# qs = [ 0 0 0 1 1 1 2 2 2 ]
		for (j in 1:length(levels)) {
		  d.subplot[i,j] <- d(communityMatrix[i,],lev=levels[j],q=qs[j])
		}
	}
	print(paste("by subplots pairs = ", rowsNum, sep=""))
	
	# between plots, not in same plot (180 pairs for 20 subplots)
	rowsNum <- nrow(communityMatrix) * (nrow(communityMatrix) + 1) / 2 - nrow(communityMatrix) - nrow(communityMatrix) / 2 		
	d.between.plot<-matrix(0,nrow=rowsNum,ncol=length(levels))
	count=0
	for(i in 1:(nrow(communityMatrix)-1)){
		for(k in (i+1):nrow(communityMatrix)){
			if ( ! (i %% 2 != 0 && k-i == 1) ) {
				count=count+1
				for (j in 1:length(levels)) {
				   d.between.plot[count,j] <- d(communityMatrix[c(i,k),],lev=levels[j],q=qs[j])
				}
			}
		}
	}
	if (count != rowsNum) stop("incorrect pairwise comparisons !")
	print(paste("between plots pairs = ", rowsNum, sep=""))

	# within plots (10 pairs)
	rowsNum <- nrow(communityMatrix) / 2
	d.plot<-matrix(0,nrow=rowsNum,ncol=length(levels))
	for(i in 1:rowsNum){	
		# d.plot[i,2] <- d(communityMatrix[c(2*i-1,2*i),],lev="gamma",q=1)
		for (j in 1:length(levels)) {
		   d.plot[i,j] <- d(communityMatrix[c(2*i-1,2*i),],lev=levels[j],q=qs[j])
		}
	}	
	print(paste("within plots pairs = ", rowsNum, sep=""))
    
    
	if (expId == 1) {	  
		
		diverBetweenPlot <- data.frame(gamma0=as.vector(d.between.plot[,1]))
		diverBetweenPlot$alpha0 <- as.vector(d.between.plot[,2])
		diverBetweenPlot$beta0 <- as.vector(d.between.plot[,3])
		diverBetweenPlot$gamma1 <- as.vector(d.between.plot[,4])
		diverBetweenPlot$alpha1 <- as.vector(d.between.plot[,5])
		diverBetweenPlot$beta1 <- as.vector(d.between.plot[,6])
		diverBetweenPlot$gamma2 <- as.vector(d.between.plot[,7])
		diverBetweenPlot$alpha2 <- as.vector(d.between.plot[,8])
		diverBetweenPlot$beta2 <- as.vector(d.between.plot[,9])
		diverBetweenPlot$region <- matrixNames[expId]
		
		diverInAPlot <- data.frame(gamma0=as.vector(d.plot[,1]))
		diverInAPlot$alpha0 <- as.vector(d.plot[,2])
		diverInAPlot$beta0 <- as.vector(d.plot[,3])
		diverInAPlot$gamma1 <- as.vector(d.plot[,4])
		diverInAPlot$alpha1 <- as.vector(d.plot[,5])
		diverInAPlot$beta1 <- as.vector(d.plot[,6])
		diverInAPlot$gamma2 <- as.vector(d.plot[,7])
		diverInAPlot$alpha2 <- as.vector(d.plot[,8])
		diverInAPlot$beta2 <- as.vector(d.plot[,9])
		diverInAPlot$region <- matrixNames[expId]
		
	} else {
		
		diverBetweenPlot_tmp <- data.frame(gamma0=as.vector(d.between.plot[,1]))
		diverBetweenPlot_tmp$alpha0 <- as.vector(d.between.plot[,2])
		diverBetweenPlot_tmp$beta0 <- as.vector(d.between.plot[,3])
		diverBetweenPlot_tmp$gamma1 <- as.vector(d.between.plot[,4])
		diverBetweenPlot_tmp$alpha1 <- as.vector(d.between.plot[,5])
		diverBetweenPlot_tmp$beta1 <- as.vector(d.between.plot[,6])
		diverBetweenPlot_tmp$gamma2 <- as.vector(d.between.plot[,7])
		diverBetweenPlot_tmp$alpha2 <- as.vector(d.between.plot[,8])
		diverBetweenPlot_tmp$beta2 <- as.vector(d.between.plot[,9])
		diverBetweenPlot_tmp$region <- matrixNames[expId]
		diverBetweenPlot <- rbind(diverBetweenPlot, diverBetweenPlot_tmp)
		
		diverInAPlot_tmp <- data.frame(gamma0=as.vector(d.plot[,1]))
		diverInAPlot_tmp$alpha0 <- as.vector(d.plot[,2])
		diverInAPlot_tmp$beta0 <- as.vector(d.plot[,3])
		diverInAPlot_tmp$gamma1 <- as.vector(d.plot[,4])
		diverInAPlot_tmp$alpha1 <- as.vector(d.plot[,5])
		diverInAPlot_tmp$beta1 <- as.vector(d.plot[,6])
		diverInAPlot_tmp$gamma2 <- as.vector(d.plot[,7])
		diverInAPlot_tmp$alpha2 <- as.vector(d.plot[,8])
		diverInAPlot_tmp$beta2 <- as.vector(d.plot[,9])	
		diverInAPlot_tmp$region <- matrixNames[expId]
		diverInAPlot <- rbind(diverInAPlot, diverInAPlot_tmp) 	
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



