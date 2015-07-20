# https://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/
library(ggplot2)
library(scales)
library(reshape2)
library(plyr)

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

n <- length(matrixNames) 
m <- length(matrixNamesNo454)

otuThr = 97

source("Modules/init.R", local=TRUE)
source("Modules/MaximizeDiversity.R", local=TRUE)

# all possible combinations of m_comb plots
m_comb = 4

######## 454 #######	

for (expId in 1:n) {	
    communityMatrix <- init(expId, otuThr, "-by-plot")
    
    print(paste(matrixNames[expId], " (", otuThr, "%) : OTUs = ", ncol(communityMatrix), ", total reads = ", sum(communityMatrix), "; min sample per plot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))
		    
	######## combinations of 4 #######
	rangeD <- getDiversitiesOfCombPlots(communityMatrix, m_comb) 
	
	maxD <- rangeD$maxD
	
#	minD <- rangeD$minD
#	print(xtable(maxD),sanitize.text.function=function(x){x})
#	print(xtable(minD),sanitize.text.function=function(x){x})

	# removed max value
	maxD <- maxD[,-ncol(maxD)]
		
	if (expId==1) {
		gamma_0 <- matrix(0,nrow=nrow(communityMatrix),ncol=(n+m))
		rownames(gamma_0) <- rownames(communityMatrix)
		colnames(gamma_0) <- c(matrixNames, matrixNamesNo454)

		beta_1 <- matrix(0,nrow=nrow(communityMatrix),ncol=(n+m))
		rownames(beta_1) <- rownames(communityMatrix)
		colnames(beta_1) <- c(matrixNames, matrixNamesNo454)
		
		rownames(gamma_0) <- gsub("CM30C30", "Plot9", rownames(gamma_0), ignore.case = T)
		rownames(gamma_0) <- gsub("LB1", "Plot10", rownames(gamma_0), ignore.case = T)
		rownames(beta_1) <- gsub("CM30C30", "Plot9", rownames(beta_1), ignore.case = T)
		rownames(beta_1) <- gsub("LB1", "Plot10", rownames(beta_1), ignore.case = T)	
	}

	if(all(tolower(rownames(gamma_0)) != tolower(colnames(maxD)))) 
		stop(paste("Plot name does not match in ", matrixNames[expId]))

	# gamma_0	
	gamma_0[,expId] <- maxD[1,]
	# beta_1	
	beta_1[,expId] <- maxD[6,]
}

######## non 454 #######	

for (expId in 1:m) {    	
    communityMatrix <- initNon454ByPlot(expId, otuThr)
        
    print(paste(matrixNamesNo454[expId], " (", otuThr, "%) : OTUs = ", ncol(communityMatrix), ", total reads = ", sum(communityMatrix), "; min sample per plot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))
		    
	######## combinations of 4 #######
	rangeD <- getDiversitiesOfCombPlots(communityMatrix, m_comb) 

	maxD <- rangeD$maxD
#	minD <- rangeD$minD
#	print(xtable(maxD),sanitize.text.function=function(x){x})
#	print(xtable(minD),sanitize.text.function=function(x){x})

	# removed max value
	maxD <- maxD[,-ncol(maxD)]
		
	if (expId==3) { # hard code for missing data
		maxD <- cbind(maxD, Plot7=rep(0,nrow(maxD)))
		maxD <- cbind(maxD, Plot8=rep(0,nrow(maxD)))
		maxD <- maxD[,c(1:6,9,10,7:8)]
	} 
	
	if(all(tolower(rownames(gamma_0)) != tolower(colnames(maxD)))) 
		stop(paste("Plot name does not match in ", matrixNamesNo454[expId]))

	# gamma_0	
	gamma_0[,(expId+n)] <- maxD[1,]
	# beta_1	
	beta_1[,(expId+n)] <- maxD[6,]
}

######## heatmap #######
divs <- c("gamma_0", "beta_1")
subTitles <- c(expression(paste("(a) "^"0", D[gamma])),expression(paste("(b) "^"1", D[beta])))

for (i in 1:length(divs)) {

	rank_table <- as.data.frame(get(divs[i]))
	#rank_table[is.na(rank_table)] <- 0

	# log
	outputRank <- paste(workingPath, "prob-4plots-", divs[i], ".txt", sep = "")
	write.table(rank_table, outputRank, sep="\t", row.names=T)

	rank_table$Plots <- rownames(rank_table)
	rank.m <- melt(rank_table)

	colnames(rank.m)[grep("^value$", colnames(rank.m))] <- "Probability" # replace "value"
	rank.m$variable <- gsub("invertebrates", "inverts", rank.m$variable, ignore.case = T)
	rank.m$Plots <- factor(rank.m$Plots, levels=rev(unique(rank.m$Plots)))
	rank.m$variable <- factor(rank.m$variable, levels=unique(rank.m$variable))

#	pdf(paste(workingPath, figDir, "/prob-4plots-", divs[i], ".pdf", sep = ""), width=6, height=6)

	mdsp <- ggplot(rank.m, aes(x=variable, y=Plots, fill = Probability)) + 
			 geom_tile(colour = "white") + 
			 scale_fill_gradient(high = "steelblue", low = "white") + guides(fill = guide_legend(reverse=TRUE)) + ggtitle(subTitles[i]) +
			 theme(legend.position="top", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) 
	
	if (i==1) mylegend<-g_legend(mdsp) 		 
 
	mdsp <- mdsp + theme(legend.position="none")
	assign(paste('mdsp', i, sep=''), mdsp)
#	invisible(dev.off())         
}

pdf(paste(workingPath, figDir, "/prob-4plots-gamma0-beta1.pdf", sep = ""), width=8, height=5)	   
grid.arrange(mylegend,arrangeGrob(mdsp1,mdsp2,ncol = 2, nrow=1), ncol=1, nrow=2, heights=c(1/15,14/15))
invisible(dev.off()) 
