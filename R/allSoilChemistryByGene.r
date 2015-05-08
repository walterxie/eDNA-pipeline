# "-by-plot" trigger merge 2 subplots columns
library(RColorBrewer)

source("init.R", local=TRUE)

beta1 <- function(communityMatrix) {    
	####### beta 1 ########
	# including diagonal
    rowsNum <- nrow(communityMatrix) * (nrow(communityMatrix) + 1) / 2  		
	d.beta1 <- matrix(0,nrow=nrow(communityMatrix),ncol=nrow(communityMatrix))
	count=0
	for(i in 1:nrow(communityMatrix)){
		for(k in i:nrow(communityMatrix)){
			count=count+1				
			d.beta1[k,i] <- d(communityMatrix[c(i,k),],lev="beta",q=1)-1				
		}
	}
	if (count != rowsNum) stop("incorrect pairwise comparisons !")
	
	return (as.dist(d.beta1))
}

sourcePath <- "~/Subversion/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
workingPath <- "~/Subversion/compevol/research/NZGenomicObservatory/Metabarcoding/"
experiments <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name
matrixNames <-  experiments
n <- length(matrixNames) 

matrixNamesNo454 <-  c("seedlings","trees","invertebrates","birds") # need expId correct for "birds","seedlings" 
m <- length(matrixNamesNo454)

otuThr = 97

inputEnvData <- paste(workingPath, "data/LJ12027.txt", sep="")    
envData <- read.table(inputEnvData, header=T, row.names=1, sep="\t")
sampleFactor <- rep(1:10, each=2)
adonis_formula <- d.beta1    

######## 454 #######

for (expId in 1:n) {	
    initList <- init(expId, otuThr, "-by-subplot")
    communityMatrix <- initList$communityMatrix

	####### beta 1 ########
	d.beta1 <- beta1(communityMatrix)
    #mantel.beta1 <- mantel(elevPlotDist, d.beta1, permutations=4999)

	aov <- adonis(formula = eval(adonis_formula) ~ pH + Olsen.P + EC + Organic.C + C.N.ratio +
		  Total.N + NO3.N + NH4.N + Water.Content, 
		  data = envData,
		  permutations = 1999, 
		  strata = as.factor(sampleFactor))
			  
    print(paste(matrixNames[expId], " (", otuThr, "%) by subplot : aov = ", sep=""))
    #print(aov$aov.tab)
    print(xtable(aov$aov.tab), type = "html")
    #print(xtable(aov$aov.tab))	
}

######## non 454 #######

for (expId in 1:m) {    	
    initList <- initNon454ByPlot(expId, otuThr)
    communityMatrix <- initList$communityMatrix
    
	####### beta 1 ########
	# (45 pairs for 10 plots)
    d.beta1 <- beta1(communityMatrix)
    #mantel.beta1 <- mantel(elevPlotDist, d.beta1, permutations=4999)
    
	envData <- read.table(inputEnvData, header=T, row.names=1, sep="\t")
	sampleFactor <- rep(1:10, each=2)
    adonis_formula <- d.beta1    
	aov <- adonis(formula = eval(adonis_formula) ~ pH + Olsen.P + EC + Organic.C + C.N.ratio +
		  Total.N + NO3.N + NH4.N + Water.Content, 
		  data = envData,
		  permutations = 1999, 
		  strata = as.factor(sampleFactor))
			  
    print(paste(matrixNamesNo454[expId], " by plot : aov = ", sep=""))
    #print(aov$aov.tab)
    print(xtable(aov$aov.tab), type = "html")
    #print(xtable(aov$aov.tab))	
}



