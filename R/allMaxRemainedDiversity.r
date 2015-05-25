
library(gplots)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(xtable)

# change config below
sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
setwd(sourcePath)

workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
experiments <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # # only for cm file name and folder name
matrixNames <-  experiments

matrixNamesNo454 <-  c("seedlings","trees","invertebrates","birds") # need expId correct for "birds","seedlings" 

n <- length(matrixNames) 
mypalette <- c("red", "orange", "green", "purple", "blue", "brown")
myshape <- seq(0, (0 + n-1))

m <- length(matrixNamesNo454)
mypalette2 <- brewer.pal(m,"Dark2")
myshape2 <- seq(0, (0 + m-1))

otuThr = 97

source("init.R", local=TRUE)
source("MaximizeDiversity.R", local=TRUE)

############### gamma 0 ##############

q=0
ylab=expression(paste("% of maximum "^"0", D[gamma]))
for (lev in c("gamma", "beta")) {
    
    if (lev == "beta") {
		q=1
		ylab=expression(paste("% of maximum "^"1", D[beta]))
	}
	

	allRanks <- NULL
	allMaxDiv <- NULL

	######## 454 #######

	for (expId in 1:n) {	
		communityMatrix <- init(expId, otuThr, "-by-plot")
	
		# 3 columns: rank, diversity, removedSites
		maxDiv <- getMaxRemainedDiversity(communityMatrix, level, q)
			
		if (expId==1) {		
			allRanks <- data.frame(Plots= maxDiv[order(maxDiv$removedSites),][,3], stringsAsFactors = FALSE)
			allMaxDiv <- data.frame(Rank= maxDiv[,1], stringsAsFactors = FALSE)		
		} else {
		   if(all(allRanks$Plots != maxDiv[order(maxDiv$removedSites),][,3])) 
			   stop(paste("Find different plots ", allRanks$Plots, " != ", maxDiv[order(maxDiv$removedSites),][,1], sep=""))
		}
	
		allRanks[,matrixNames[expId]] <- as.numeric(maxDiv[order(maxDiv$removedSites),][,1]) 
		allMaxDiv[,matrixNames[expId]] <- as.numeric(maxDiv[,2])
	}

	######## non 454 #######	

	for (expId in 1:m) {    	
		communityMatrix <- initNon454ByPlot(expId, otuThr)

		rownames(communityMatrix) <- gsub("CM30C30", "Plot9", rownames(communityMatrix), ignore.case = T)
		rownames(communityMatrix) <- gsub("LB1", "Plot10", rownames(communityMatrix), ignore.case = T)

		# 3 columns: rank, diversity, removedSites
		maxDiv <- getMaxRemainedDiversity(communityMatrix, level, q)
	
		#"invertebrates", missing Plot7, Plot8
		if (expId==3) {
		   newrow <- c(NA, 0, "Plot7")
		   maxDiv <- rbind(maxDiv,newrow)
		   newrow <- c(NA, 0, "Plot8")
		   maxDiv <- rbind(maxDiv,newrow)
		}

		if(all(allRanks$Plots != maxDiv[order(maxDiv$removedSites),][,3])) 
			stop(paste("Find different plots ", allRanks$Plots, " != ", maxDiv[order(maxDiv$removedSites),][,1], sep=""))

		allRanks[,matrixNamesNo454[expId]] <- as.numeric(maxDiv[order(maxDiv$removedSites),][,1])
		allMaxDiv[,matrixNamesNo454[expId]] <- as.numeric(maxDiv[,2])
	}

	######## file #######
	outputMaxDiv <- paste(workingPath, "plots-rank-table-", lev, q, ".txt", sep = "")
	write.table(allRanks, outputMaxDiv, sep="\t", row.names=FALSE)
	#write.csv(allRanks, outputMaxDiv, row.names=FALSE)

	######## correlation #######
	corSp <- cor(allRanks[,-1], use ="pairwise.complete.obs", method="spearman")
	# set upper tri to NA
	corSp[upper.tri(corSp)] <- NA
	corSp[corSp=="1"] <- NA
	corSp <- corSp[-1,-ncol(corSp)]
	# take significate digits 2
	corSp <- formatC(signif(corSp,digits=2), digits=2,format="fg", flag="#")
	
	outputCor <- paste(workingPath, "plots-rank-correlation-", lev, q, ".txt", sep = "")
	write.table(corSp, outputCor, na="", sep="\t")
#	print(xtable(corSp))
#	cor(allMaxDiv[,-which(names(allMaxDiv) == "Rank" | names(allMaxDiv) == "invertebrates")], use ="everything", method="spearman")
	
	######## correlation significance #######
	sigSp <- matrix(0,nrow=(ncol(allRanks)-1),ncol=(ncol(allRanks)-1)) 
	colnames(sigSp) <- colnames(allRanks)[-1]
	rownames(sigSp) <- colnames(allRanks)[-1]
	for (i in 2:(ncol(allRanks)-1)) { # 1st col is plots name
		for (j in (i+1):ncol(allRanks)) {	     
			tmpDf<-allRanks[,c(i,j)]
			tmpDf<-tmpDf[complete.cases(tmpDf),] # remove NA rows pairwised
			sigSp[j-1,i-1] <- cor.test(tmpDf[,1], tmpDf[,2], method="spearman")$p.value
		}
	}
	
	sigSp[upper.tri(sigSp)] <- NA
	sigSp[sigSp=="0"] <- NA
	sigSp <- sigSp[-1,-ncol(sigSp)]
	sigSp <- formatC(signif(sigSp,digits=2), digits=2,format="fg", flag="#")
	
	outputSig <- paste(workingPath, "plots-rank-cor-sig-", lev, q, ".txt", sep = "")
	write.table(sigSp, outputSig, na="", sep="\t")
	
	######## all in 1 table #######
	allSp <- matrix(0,nrow=nrow(corSp),ncol=ncol(corSp)) 
	colnames(allSp) <- colnames(corSp)
	rownames(allSp) <- rownames(corSp)

	for (i in 1:ncol(allSp)) {
		for (j in 1:nrow(allSp)) {	     
			allSp[i,j] <- paste(corSp[i,j], " (", sigSp[i,j], ")", sep = "")
		}
	}
	
	print(xtable(allSp, caption = "Spearman correlations and significance", label = paste("tab:corMaxRemained", lev, q, sep = "")), 
			sanitize.text.function = function(x){x})

	
	######## figure #######
    allMaxDiv <- allMaxDiv[,c(2:ncol(allMaxDiv),1)] # move Rank col to the last col

	pdf(paste(workingPath, "figures/max-remained-", lev, q, ".pdf", sep = ""), width=9, height=5)
	attach(mtcars)
	par(mfrow=c(1,2), oma=c(4,0,0,0)) 

	# 1st col is Rank
	for (expId in 1:n) {
	    yVal = allMaxDiv[,expId]/max(allMaxDiv[,expId])	
		if (expId == 1) {
			par(mar=c(4,5,2,1))		
			plot(allMaxDiv$Rank, yVal, pch=myshape[expId], col=mypalette[expId], xlim=c(nrow(allMaxDiv),1), ylim=c(0,1), 
			   xlab="number of sites remained", ylab=ylab, main="(a)") 				  		
		} else {
			points(allMaxDiv$Rank, yVal, pch=myshape[expId], col=mypalette[expId]) 
		}

		lines(allMaxDiv$Rank, yVal, lty=2, col=mypalette[expId]) 	
	}

	for (expId in 1:m) {
		incr = n
		yVal = allMaxDiv[,expId+incr]/max(allMaxDiv[,expId+incr])
		if (expId == 1) {
			par(mar=c(4,5,2,1))		
			plot(allMaxDiv$Rank, yVal, pch=myshape2[expId], col=mypalette2[expId], xlim=c(nrow(allMaxDiv),1), ylim=c(0,1), 
			   xlab="number of sites remained", ylab=ylab, main="(b)") 				  		
		} else if (expId == 3) { # 2 sites missing
			points(allMaxDiv$Rank[-c(1,2)], yVal[-c(9,10)], pch=myshape2[expId], col=mypalette2[expId]) 
			lines(allMaxDiv$Rank[-c(1,2)], yVal[-c(9,10)], lty=2, col=mypalette2[expId]) 
		} else {
		    points(allMaxDiv$Rank, yVal, pch=myshape2[expId], col=mypalette2[expId]) 
		}
		
		if (expId != 3) {
			lines(allMaxDiv$Rank, yVal, lty=2, col=mypalette2[expId]) 
		}
	}

	par(usr=c(0,1,0,1),xpd=NA)  
	legend(-1.38, -0.3, ncol=3, legend=matrixNames, pch=as.numeric(myshape), col=mypalette)
	legend(0.1, -0.3, ncol=2, legend=matrixNamesNo454, pch=as.numeric(myshape2), col=mypalette2)
	invisible(dev.off())  
}


