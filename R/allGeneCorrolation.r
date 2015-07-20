library(vegan)
library(vegetarian)
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

n <- length(matrixNames) 
m <- length(matrixNamesNo454)

otuThr = 97

source("Modules/init.R", local=TRUE)

plots <- c("Plot1","Plot2","Plot3","Plot4","Plot5","Plot6","Plot7","Plot8","Plot9","Plot10")

######## 454 #######
beta1DisList <- list()

for (expId in 1:n) {	
    communityMatrix <- init(expId, otuThr, "-by-plot")
	
	if ( all( rownames(communityMatrix) != plots ) )
		stop( paste("community matrix", matrixNames[expId], "plot names are incorrect !") )

	####### beta 1 ########
	# (45 pairs for 10 plots)
    d.beta1 <- beta1(communityMatrix)
    
    beta1DisList[[ expId ]] <- d.beta1
	    
    print(paste(matrixNames[expId], " (", otuThr, "%) : total = ", sum(communityMatrix), "; min sample per subplot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))
	
}

######## non 454 #######	
beta1No454DisList <- list()

for (expId in 1:m) {    	
    communityMatrix <- initNon454ByPlot(expId, otuThr)
	rownames(communityMatrix) <- gsub("CM30C30", "Plot9", rownames(communityMatrix), ignore.case = T)
	rownames(communityMatrix) <- gsub("LB1", "Plot10", rownames(communityMatrix), ignore.case = T)

	if (matrixNamesNo454[expId] == "invertebrates") {
		if ( all( rownames(communityMatrix) != plots[-c(7,8)] ) )
			stop( paste("community matrix", matrixNamesNo454[expId], "plot names are incorrect !") )
	} else {
		if ( all( rownames(communityMatrix) != plots ) )
			stop( paste("community matrix", matrixNamesNo454[expId], "plot names are incorrect !") )
	}

	####### beta 1 ########
	# (45 pairs for 10 plots)
    d.beta1 <- beta1(communityMatrix)
    
    beta1No454DisList[[ expId ]] <- d.beta1
       
    print(paste(matrixNamesNo454[expId], " (", otuThr, "%) : total = ", sum(communityMatrix), "; min sample per subplot = ", min(rowSums(communityMatrix)), "; max = ", max(rowSums(communityMatrix)), sep=""))

}

m.mantel <- matrix(0,nrow=(n+m),ncol=(n+m))
colnames(m.mantel) <- c(matrixNames, matrixNamesNo454)
rownames(m.mantel) <- c(matrixNames, matrixNamesNo454)
m.signif <- matrix(0,nrow=(n+m),ncol=(n+m))
colnames(m.signif) <- c(matrixNames, matrixNamesNo454)
rownames(m.signif) <- c(matrixNames, matrixNamesNo454)

for (expId1 in 1:(n+m-1)) {
	for (expId2 in (expId1+1):(n+m)) {	     
		if (expId1 <= n) {
			d.comm1 <- beta1DisList[[ expId1 ]]
			exp.comm1 <- matrixNames[expId1]
		} else {
			d.comm1 <- beta1No454DisList[[ expId1-n ]]
			exp.comm1 <- matrixNamesNo454[expId1-n]
		}

		if (expId2 <= n) {
			d.comm2 <- beta1DisList[[ expId2 ]]
			exp.comm2 <- matrixNames[expId2]
		} else {
			d.comm2 <- beta1No454DisList[[ expId2-n ]]
			exp.comm2 <- matrixNamesNo454[expId2-n]
		}

		# hard code for missing data in invertebrates
		if (length(d.comm1) > length(d.comm2)) {
			m.comm <- as.matrix(d.comm1) 
			d.comm1 <- as.dist(m.comm[c(-7,-8),c(-7,-8)]) 
		} else if (length(d.comm1) < length(d.comm2)) {
			m.comm <- as.matrix(d.comm2) 
			d.comm2 <- as.dist(m.comm[c(-7,-8),c(-7,-8)]) 	     
		}

		####### mantel test #######
		mantel.comm <- mantel(d.comm1, d.comm2, permutations=4999) 

		m.mantel[expId2, expId1] <- mantel.comm$statistic
		m.mantel[expId1, expId2] <- mantel.comm$statistic
		m.signif[expId2, expId1] <- mantel.comm$signif
		m.signif[expId1, expId2] <- mantel.comm$signif

		print(paste(exp.comm1, " vs. ", exp.comm2, ", mantel statistic r = ", round(mantel.comm$statistic, 3), "; significance = ", mantel.comm$signif, "; permutations = ", mantel.comm$permutations, sep=""))
	}
}

####### MDS of correlation #####
# Classical MDS
for (expId1 in 1:(n+m)) {
    m.mantel[expId1, expId1] <- 1
}
d.mantel <- 1-m.mantel # euclidean distances between the rows

fit <- cmdscale(d.mantel, eig=TRUE, k=2) # k is the number of dim

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]

labels = gsub("invertebrates", "inverts", row.names(fit$points), ignore.case = T)

pdf(paste(workingPath, figDir, "/mds-pairewise-cm-", otuThr, ".pdf", sep = ""), width=5, height=5)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS (1-Correlation)",	type="n")
text(x, y, labels = labels, cex=.7, xpd=TRUE)
invisible(dev.off()) 

output <- paste(workingPath, "mds-pairewise-mantel-", otuThr, ".csv", sep = "")
# tsv cannot display the 1st cell of columns
write.csv(d.mantel, output, quote=FALSE)

# no birds
m.mantel1 <- m.mantel[-(n+m),-(n+m)] 
d.mantel <- 1-m.mantel1 # euclidean distances between the rows
fit <- cmdscale(d.mantel,eig=TRUE, k=2) # k is the number of dim

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]

labels = gsub("invertebrates", "inverts", row.names(fit$points), ignore.case = T)

pdf(paste(workingPath, figDir, "/mds-no-birds-pairewise-cm-", otuThr, ".pdf", sep = ""), width=5, height=5)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS (1-Correlation)", type="n")
text(x, y, labels = labels, cex=.7, xpd=TRUE)
invisible(dev.off()) 

####### mantel test #######
df.mantel <- data.frame(m.mantel, check.names=FALSE)
for (expId1 in 1:(n+m-1)) {
	for (expId2 in (expId1+1):(n+m)) {	     
        df.mantel[expId1, expId2] <- c("")
        tmp <- paste(round(m.mantel[expId2, expId1], 3), " (", m.signif[expId2, expId1], ")", sep="")
        df.mantel[expId2, expId1] <- c(tmp)
	}
}
for (expId1 in 1:(n+m)) {
    df.mantel[expId1, expId1] <- c("")
}
df.mantel <- df.mantel[-1,-ncol(df.mantel)]
 
print(xtable(df.mantel, caption = "Pairwise community matrix correlations of effective $\\beta$ diversity within and between 
	the eDNA  data sets and traditional data sets, Mantel statistic $r$, and their significance in parentheses 
	using Mantel's test based on 4999 permutations.", label = "tab:geneCorr", caption.placement = "top"), 
	sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

