
library(vegan)
library(vegetarian)
library(data.table)
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
if(!exists("levels")) stop("levels of Jost diversity are missing !")
if(!exists("qs")) stop("qs of Jost diversity are missing !")

n <- length(matrixNames) 
m <- length(matrixNamesNo454)

allMatrixNames <-  c(matrixNames, matrixNamesNo454)

otuThr = 97

source("Modules/init.R", local=TRUE)

plots <- c("Plot1","Plot2","Plot3","Plot4","Plot5","Plot6","Plot7","Plot8","Plot9","Plot10")

cmdsList <- list()
cmds8PlotsList <- list()

######## 454 and non 454 #######	

for (expId in 1:(n+m)) {
	if (expId > n) {	
		communityMatrix <- initNon454ByPlot(expId-n, otuThr)
    } else {
		communityMatrix <- init(expId, otuThr, "-by-plot")  
    }
	rownames(communityMatrix) <- gsub("CM30C30", "Plot9", rownames(communityMatrix), ignore.case = T)
	rownames(communityMatrix) <- gsub("LB1", "Plot10", rownames(communityMatrix), ignore.case = T)
	
	if (allMatrixNames[expId] == "invertebrates") {
		if ( all( rownames(communityMatrix) != plots[-c(7,8)] ) )
			stop( paste("community matrix", allMatrixNames[expId], "plot names are incorrect !") )
	} else {
		if ( all( rownames(communityMatrix) != plots ) )
			stop( paste("community matrix", allMatrixNames[expId], "plot names are incorrect !") )
	}
	
    d.beta1 <- beta1(communityMatrix)
		
	# Classical multidimensional scaling
	cmds <- cmdscale(d.beta1, eig=TRUE, k=2) # k is the number of dim
	print(cmds)

	#### 8 plots mds for Procrustes
	if (allMatrixNames[expId] == "invertebrates") {
		cmds8PlotsList[[ expId ]] <- cmds	
	} else {
		cmdsList[[ expId ]] <- cmds
	
		d.beta1.8p <- beta1(communityMatrix[-c(7,8),]) # hard code for missing data
		cmds8Plots <- cmdscale(d.beta1.8p, eig=TRUE, k=2)
		cmds8PlotsList[[ expId ]] <- cmds8Plots	
	}	
		   
} #END for expId

######## Procrustes analysis #######
proc <- matrix(0,nrow=(n+m),ncol=(n+m))
colnames(proc) <- allMatrixNames
rownames(proc) <- allMatrixNames

prot <- matrix(0,nrow=(n+m),ncol=(n+m))
colnames(prot) <- allMatrixNames
rownames(prot) <- allMatrixNames
prot.signif <- matrix(0,nrow=(n+m),ncol=(n+m))
colnames(prot.signif) <- allMatrixNames
rownames(prot.signif) <- allMatrixNames

for (i in 1:(n+m-1)) {
	for (j in (i+1):(n+m)) {	
		if (allMatrixNames[i] == "invertebrates" || allMatrixNames[j] == "invertebrates") {
			cmds1 <- cmds8PlotsList[[ i ]]
			cmds2 <- cmds8PlotsList[[ j ]]
		} else {
			cmds1 <- cmdsList[[ i ]]
			cmds2 <- cmdsList[[ j ]]
		}

		# symmetric Procrustes plot.
		vare.proc <- procrustes(cmds1, cmds2)
		
		proc[j, i] <- vare.proc$ss # Sum of Squares: lower triangle

		vare.proc <- procrustes(cmds2, cmds1)
		proc[i, j] <- vare.proc$ss

		# repeatedly to estimate the "significance" of the Procrustes statistic.
		p<-protest(cmds1, cmds2, scores = "sites", permutations = 4999)
		
		prot[j, i] <- p$ss # Sum of Squares: lower triangle
		prot[i, j] <- p$scale # Correlation in a symmetric Procrustes rotation: upper triangle
				
		prot.signif[j, i] <- p$signif # lower triangle
				
	} #END for j
} #END for i


######## Procrustes plot #######

subT=1
pdf(paste(workingPath, figDir, "/procrustes-all.pdf", sep = ""), width=7, height=10)
attach(mtcars)
par(mfrow=c(5,3), oma=c(0,1,0,0), mai = c(0.3,0.2,0.2,0.2))

for (i in 1:(n+m-1)) {
	for (j in (i+1):(n+m)) {	
		if (allMatrixNames[i] == "invertebrates" || allMatrixNames[j] == "invertebrates") {
			cmds1 <- cmds8PlotsList[[ i ]]
			cmds2 <- cmds8PlotsList[[ j ]]
		} else {
			cmds1 <- cmdsList[[ i ]]
			cmds2 <- cmdsList[[ j ]]
		}

		# symmetric Procrustes plot.
		vare.proc <- procrustes(cmds1, cmds2, symmetric = TRUE)
		labels = gsub("Plot", "", rownames(cmds1$points), ignore.case = T)
				
		if (i < 6 && (j==7 || j==8 || j==9)) {
			sig <- ""
			if (prot.signif[j, i] < 0.005) {
				sig <- "(***)"
			} else if (prot.signif[j, i] < 0.01) {
				sig <- "(**)"
			} else if (prot.signif[j, i] < 0.05) {
				sig <- "(*)"
			}
			main=bquote("("*.(letters[subT])*")"~.(allMatrixNames[i])%<-%.(allMatrixNames[j])~.(sig))
#			main=bquote("("*.(letters[subT])*")"~.(allMatrixNames[i])%<-%.(allMatrixNames[j])~" (ss="~.(round(vare.proc$ss, 2))~")")
			subT=subT+1
			plot(vare.proc, main=main, ylab="", xlab="")
			text(vare.proc, labels=labels, cex=1, xpd=TRUE)						
		} 
						
	} #END for j
} #END for i
invisible(dev.off())

######## Procrustes statistic #######
prot.ss <- prot
prot.ss[upper.tri(prot.ss)] <- 0
prot.ss <- formatC(signif(prot.ss,digits=2), digits=2,format="fg", flag="#")
prot.signif[upper.tri(prot.signif)] <- 0
prot.signif <- formatC(signif(prot.signif,digits=2), digits=2,format="fg", flag="#")

prot.table <- matrix( paste(prot.ss, " (", prot.signif, ")", sep=""), nrow=nrow(prot.ss), dimnames=dimnames(prot.ss) )
prot.table[prot.table=="0 (0)"] <- ""

prot.table <- prot.table[-1,-ncol(prot.table)]

print(xtable(prot.table, caption = "Pairwise community matrix Procrustes analysis of effective $\\beta$ diversity within and between 
	the eDNA data sets and traditional data sets, estimated sum of squared differences and their significance in parentheses 
	based on 4999 permutations.", label = "tab:protest", caption.placement = "top"), 
	sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

######## NMDS plot of ss #######
proc[upper.tri(proc)] <- 0

fit <- cmdscale(proc, eig=TRUE, k=2) # k is the number of dim

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]

labels = gsub("invertebrates", "inverts", row.names(fit$points), ignore.case = T)

pdf(paste(workingPath, figDir, "/mds-pairewise-procrustes-", otuThr, ".pdf", sep = ""), width=5, height=5)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS (Sum of Squared Differences)", type="n")
text(x, y, labels = labels, cex=.7, xpd=TRUE)
invisible(dev.off()) 

output <- paste(workingPath, "mds-pairewise-procrustes-", otuThr, ".csv", sep = "")
# tsv cannot display the 1st cell of columns
write.csv(proc, output, quote=FALSE)

proc1 <- proc[-(n+m),-(n+m)] # no birds
fit <- cmdscale(proc1, eig=TRUE, k=2) # k is the number of dim

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]

labels = gsub("invertebrates", "inverts", row.names(fit$points), ignore.case = T)

pdf(paste(workingPath, figDir, "/mds-no-birds-pairewise-procrustes-", otuThr, ".pdf", sep = ""), width=5, height=5)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS (Sum of Squared Differences)", type="n")
text(x, y, labels = labels, cex=.7, xpd=TRUE)
invisible(dev.off()) 

