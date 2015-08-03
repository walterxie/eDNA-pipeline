
library(gplots)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(xtable)

# change config below
#figDir <- figDir, ""
#sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
#setwd(sourcePath)
#workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
#matrixNames <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name   
#matrixNamesNo454 <-  c("seedlings","trees","invertebrates","birds") # need expId correct for "birds","seedlings" 

if(!exists("tableFile")) stop("table file is missing !")
if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")
if(!exists("matrixNamesNo454")) stop("matrix names of traditional methods are missing !")
if(!exists("rmSingleton")) stop("rmSingleton flag is missing !")

n <- length(matrixNames) 
mypalette <- c("red", "orange", "green", "purple", "blue", "brown")
myshape <- seq(0, (0 + n-1))

m <- length(matrixNamesNo454)
mypalette2 <- brewer.pal(m,"Dark2")
myshape2 <- seq(0, (0 + m-1))

otuThr = 97

source("Modules/init.R", local=TRUE)
source("Modules/MaximizeDiversity.R", local=TRUE)

############### gamma 0 ##############

q=0
ylab=expression(paste("% of maximum "^"0", D[gamma]))
for (lev in c("gamma", "beta")) {
    
    if (lev == "beta") {
		q=1
		ylab=expression(paste("% of maximum "^"1", D[beta]))
	}
	
	# ranks
	allRanks <- NULL 
	# max diversity
	allMaxDiv <- NULL

	######## 454 #######
	for (expId in 1:n) {	
		communityMatrix <- init(expId, otuThr, "-by-plot")
	
		# 3 columns: rank, diversity, site
		maxDiv <- getMaxRemainedDiversity(communityMatrix, level, q)
			
		if (expId==1) {		
			allRanks <- data.frame(row.names=rownames(maxDiv), stringsAsFactors = FALSE)
			allMaxDiv <- data.frame(row.names=rownames(maxDiv), stringsAsFactors = FALSE)		
		} else {
		   if(!all(rownames(allRanks) == rownames(maxDiv))) 
			   stop(paste("Find different plots ", rownames(allRanks), " != ", rownames(maxDiv), sep=""))
		}
	
		allRanks[,matrixNames[expId]] <- as.numeric(maxDiv[,1]) 
		allMaxDiv[,matrixNames[expId]] <- as.numeric(maxDiv[,2])
	}  # END for expId

	######## non 454 #######	
	for (expId in 1:m) {    	
		communityMatrix <- initNon454ByPlot(expId, otuThr)

		rownames(communityMatrix) <- gsub("CM30C30", "Plot9", rownames(communityMatrix), ignore.case = T)
		rownames(communityMatrix) <- gsub("LB1", "Plot10", rownames(communityMatrix), ignore.case = T)

		# 3 columns: rank, diversity, site
		maxDiv <- getMaxRemainedDiversity(communityMatrix, level, q)
	
		#"invertebrates", missing Plot7, Plot8
		if (expId==3) {
		   newrow <- c(NA, 0)
		   maxDiv <- rbind(maxDiv,newrow)
		   rownames(maxDiv)[nrow(maxDiv)] <- "Plot7"
		   newrow <- c(NA, 0)
		   maxDiv <- rbind(maxDiv,newrow)
		   rownames(maxDiv)[nrow(maxDiv)] <- "Plot8"
		   maxDiv <- maxDiv[order(rownames(maxDiv)),]
		}

		if(!all(rownames(allRanks) == rownames(maxDiv))) 
			stop(paste("Find different plots ", rownames(allRanks), " != ", rownames(maxDiv), sep=""))

		allRanks[,matrixNamesNo454[expId]] <- as.numeric(maxDiv[,1])
		allMaxDiv[,matrixNamesNo454[expId]] <- as.numeric(maxDiv[,2])
	}  # END for expId

	######## rank table #######
	allRanks <- allRanks[c(1,3:nrow(allRanks),2),]

	if (lev == "gamma") {
		print(xtable(allRanks, caption = "Ranking sampling plots by removing plots sequentially so as to 
				minimize the loss of overall $\\gamma$ diversity among the remaining plots. 
				1 is the most important and removed at the last, 10 is the least important and removed in the beginning.", 
				label = "tab:maxRemainedGamma0", caption.placement = "top"), 
				sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)
	} else {
		print(xtable(allRanks, caption = "Ranking sampling plots by removing plots sequentially so as to 
				maximize the resulting effective $\\beta$ diversity among the remaining plots. 
				1 is the most important and removed at the last, 10 is the least important and removed in the beginning.",
				label = "tab:maxRemainedBeta1", caption.placement = "top"), 
				sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)		
	}

	# average table
	avg_table <- data.frame(row.names=rownames(allRanks))

	avg_table$eDNA <- rowMeans(allRanks[,1:6], na.rm=T)
	avg_table$Trad <- rowMeans(allRanks[,7:10], na.rm=T)
	avg_table$Trad.No.Birds <- rowMeans(allRanks[,7:9], na.rm=T)
	avg_table$All <- rowMeans(allRanks, na.rm=T)
	avg_table$All.No.Birds <- rowMeans(allRanks[,1:9], na.rm=T)

	avg_table <- formatC(signif(as.matrix(avg_table),digits=2), digits=2,format="fg", flag="#")

	# std table
	std_table <- matrix(0,nrow=nrow(avg_table),ncol=ncol(avg_table))
	colnames(std_table) <- colnames(avg_table)
	rownames(std_table) <- rownames(avg_table)

	for (r in 1:nrow(std_table)) {
		std_table[r,1] <- sd(allRanks[r,1:6], na.rm=T)
		std_table[r,2] <- sd(allRanks[r,7:10], na.rm=T)
		std_table[r,3] <- sd(allRanks[r,7:9], na.rm=T)
		std_table[r,4] <- sd(allRanks[r,], na.rm=T)
		std_table[r,5] <- sd(allRanks[r,1:9], na.rm=T)
	}

	std_table <- formatC(signif(std_table,digits=3), digits=3,format="fg", flag="#")

	# final
	final_table <- matrix(0,nrow=nrow(avg_table),ncol=ncol(avg_table))
	colnames(final_table) <- colnames(avg_table)
	rownames(final_table) <- rownames(avg_table)
	for (r in 1:nrow(avg_table)) {
		final_table[r,] <- paste(avg_table[r,], " $\\pm$ ",  std_table[r,])
	}

	if (lev == "gamma") {
		print(xtable(final_table, caption = "Means and standard deviations of Table~\\ref{tab:maxRemainedGamma0} between eDNA data sets 
				for ranking plots to minimize the loss of overall $\\gamma$ diversity among the remaining plots.",
				label = "tab:avgMaxRemainedGamma0", caption.placement = "top"), 
				sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)
	} else {
		print(xtable(final_table, caption = "Means and standard deviations of Table~\\ref{tab:maxRemainedBeta1} between eDNA data sets 
				for ranking plots to maximize effective $\\beta$ diversity among the remaining plots.",
				label = "tab:avgMaxRemainedBeta1", caption.placement = "top"), 
				sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)
	}

	######## correlation and significance #######
	corSp <- cor(allRanks, use ="pairwise.complete.obs", method="spearman")
	# set upper tri to NA
	corSp[upper.tri(corSp)] <- NA
	corSp[corSp=="1"] <- NA
	corSp <- corSp[-1,-ncol(corSp)]
	# take significate digits 2
	corSp <- formatC(signif(corSp,digits=2), digits=2,format="fg", flag="#")
		
	sigSp <- matrix(0,nrow=ncol(allRanks),ncol=ncol(allRanks)) 
	colnames(sigSp) <- colnames(allRanks)
	rownames(sigSp) <- colnames(allRanks)
	for (i in 1:(ncol(allRanks)-1)) { # 1st col is plots name
		for (j in (i+1):ncol(allRanks)) {	     
			tmpDf<-allRanks[,c(i,j)]
			tmpDf<-tmpDf[complete.cases(tmpDf),] # remove NA rows pairwised
			sigSp[j,i] <- cor.test(tmpDf[,1], tmpDf[,2], method="spearman")$p.value
		}
	}
	
	sigSp[upper.tri(sigSp)] <- NA
	sigSp[sigSp=="0"] <- NA
	sigSp <- sigSp[-1,-ncol(sigSp)]
	sigSp <- formatC(signif(sigSp,digits=2), digits=2,format="fg", flag="#")
	
	######## all in 1 table #######
	if(!all(rownames(corSp) == rownames(sigSp)) || !all(colnames(corSp) == colnames(sigSp))) 
		stop(paste("Correlation matrix not match significance matrix ! "))

	allSp <- matrix(0,nrow=nrow(corSp),ncol=ncol(corSp)) 
	colnames(allSp) <- colnames(corSp)
	rownames(allSp) <- rownames(corSp)

	for (i in 1:ncol(allSp)) {
		for (j in 1:nrow(allSp)) {	     
			allSp[i,j] <- paste(corSp[i,j], " (", sigSp[i,j], ")", sep = "")
		}
	}
	
	allSp[allSp==" NA ( NA)"] <- ""
	
	if (lev == "gamma") {
		print(xtable(allSp, caption = "Spearman correlations and their significance in parentheses of 
				Table~\\ref{tab:maxRemainedGamma0} between eDNA data sets 
				for ranking plots to minimize the loss of overall $\\gamma$ diversity among the remaining plots.",
				label = "tab:corMaxRemainedGamma0", caption.placement = "top"), 
				sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)
	} else {
		print(xtable(allSp, caption = "Spearman correlations and their significance in parentheses of 
				Table~\\ref{tab:maxRemainedBeta1} between eDNA data sets 
				for ranking plots to maximize effective $\\beta$ diversity among the remaining plots.",
				label = "tab:corMaxRemainedBeta1", caption.placement = "top"), 
				sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)
	}

	######## figure #######
	allMaxDiv <- allMaxDiv[c(1,3:nrow(allMaxDiv),2),]
	
	if(!all(rownames(allRanks) == rownames(allMaxDiv))) 
		stop(paste("Find different plots ", rownames(allRanks), " != ", rownames(allMaxDiv), sep=""))
			
	pdf(paste(workingPath, figDir, "/max-remained-", lev, q, ".pdf", sep = ""), width=9, height=5)
	attach(mtcars)
	par(mfrow=c(1,2), oma=c(4,0,0,0)) 

	# 1st col is Rank
	for (expId in 1:n) {
		val <- data.frame(x=allRanks[,expId], y=allMaxDiv[,expId]/max(allMaxDiv[,expId]), stringsAsFactors = FALSE)
	    val <- val[order(val$x, decreasing=T),]
		if (expId == 1) {
			par(mar=c(4,5,2,1))		
			plot(val, pch=myshape[expId], col=mypalette[expId], xlim=c(nrow(val),1), ylim=c(0,1), 
			   xlab="number of sites remained", ylab=ylab, main="(a)") 				  		
		} else {
			points(val, pch=myshape[expId], col=mypalette[expId]) 
		}

		lines(val, lty=2, col=mypalette[expId]) 	
	}  # END for expId

	for (expId in 1:m) {
		incr = n
		val <- data.frame(x=allRanks[,expId+incr], y=allMaxDiv[,expId+incr]/max(allMaxDiv[,expId+incr]), stringsAsFactors = FALSE)
		val <- val[complete.cases(val),]
	    val <- val[order(val$x, decreasing=T),]
		if (expId == 1) {
			par(mar=c(4,5,2,1))		
			plot(val, pch=myshape2[expId], col=mypalette2[expId], xlim=c(nrow(val),1), ylim=c(0,1), 
			   xlab="number of sites remained", ylab=ylab, main="(b)") 				  		
		} else {
		    points(val, pch=myshape2[expId], col=mypalette2[expId]) 
		}
		
		lines(val, lty=2, col=mypalette2[expId]) 
	} # END for expId

	par(usr=c(0,1,0,1),xpd=NA)  
	legend(-1.38, -0.3, ncol=3, legend=matrixNames, pch=as.numeric(myshape), col=mypalette)
	legend(0.1, -0.3, ncol=2, legend=matrixNamesNo454, pch=as.numeric(myshape2), col=mypalette2)
	invisible(dev.off())  
} # END for lev


