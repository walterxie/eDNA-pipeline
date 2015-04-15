# "-by-plot" trigger merge 2 subplots columns
library(ggplot2)
library(vegan)
library(vegetarian)
library(grid)
library(gridExtra)

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

# change config below
sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
setwd(sourcePath)

workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
experiments <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name
matrixNames <-  experiments
subTitles <- c("(a)","(b)","(c)","(d)","(e)","(f)")

n <- length(matrixNames) 
mypalette <- c("red", "orange", "green", "purple", "blue", "brown")
myshape <- seq(0, (0 + n-1))

otuThr = 97

source("init.R", local=TRUE)
  
#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


######## beta 1 #######

for (expId in 1:n) {	
    # "-by-plot" trigger merge 2 subplots columns
    communityMatrix <- init(expId, otuThr, "-by-subplot")    
    matrixName <- matrixNames[expId]
    d.beta1 <- beta1(communityMatrix)
		
    # Run metaMDS, get points and stress
	mds <- metaMDS(d.beta1)
	pts_mds <- as.data.frame(mds$points)
	stress_mds <- mds$stress
	
	# Get sample labels/factors
	pts_mds$plots <- sapply(strsplit(rownames(pts_mds), "-"), "[", 1)
	# convert plot names
	pts_mds$plots <- gsub("CM30C30", "Plot9", pts_mds$plots)
	pts_mds$plots <- gsub("LB1", "Plot10", pts_mds$plots)
	rownames(pts_mds) <- gsub("CM30C30", "Plot9", rownames(pts_mds))
	rownames(pts_mds) <- gsub("LB1", "Plot10", rownames(pts_mds))

	subTitle <- paste(subTitles[expId], " ", matrixNames[expId], " (stress ", round(stress_mds, 2),")", sep = "")
	
	pts_mds$plots <- factor(pts_mds$plots,levels = unique(pts_mds$plots))
	# Plot MDS ordination
	mdsp <- ggplot(pts_mds, aes(x = MDS1, y = MDS2, color = plots)) + 
			geom_point(size = 3) + 
			geom_text(aes(label = rownames(pts_mds)), size = 3, vjust = 2) +
			geom_polygon(aes(group = pts_mds$plots), fill = NA) +
			theme(legend.position="top", legend.title=element_blank()) +
			ggtitle(subTitle) 

   assign(paste('mdsp', expId, sep=''), mdsp)
}

mdsp1 <- mdsp1 + theme(axis.title.x=element_blank(), legend.position="none")
mdsp2 <- mdsp2 + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none")
mdsp3 <- mdsp3 + theme(axis.title.x=element_blank(), legend.position="none")
mdsp4 <- mdsp4 + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none")
mdsp5 <- mdsp5 + theme(legend.position="none")
mdsp6 <- mdsp6 + theme(axis.title.y=element_blank(), legend.position="none")

pdf(paste(workingPath, "figures/mds-subplots-beta1-", otuThr, ".pdf", sep = ""), width=8, height=12)	   
mylegend<-g_legend(mdsp)
grid.arrange(mylegend,arrangeGrob(mdsp1,mdsp2,mdsp3,mdsp4,mdsp5,mdsp6,ncol = 2, nrow=3), ncol=1, nrow=2, heights=c(1/20,19/20))
invisible(dev.off()) 


######## Bray-Curtis #######

for (expId in 1:n) {	
    # "-by-plot" trigger merge 2 subplots columns
    communityMatrix <- init(expId, otuThr, "-by-subplot")    
    matrixName <- matrixNames[expId]
    d.brayBin <- vegdist(communityMatrix, method="bray", binary=TRUE)
		
    # Run metaMDS, get points and stress
	mds <- metaMDS(d.brayBin)
	pts_mds <- as.data.frame(mds$points)
	stress_mds <- mds$stress
	
	# Get sample labels/factors
	pts_mds$plots <- sapply(strsplit(rownames(pts_mds), "-"), "[", 1)
	# convert plot names
	pts_mds$plots <- gsub("CM30C30", "Plot9", pts_mds$plots)
	pts_mds$plots <- gsub("LB1", "Plot10", pts_mds$plots)
	rownames(pts_mds) <- gsub("CM30C30", "Plot9", rownames(pts_mds))
	rownames(pts_mds) <- gsub("LB1", "Plot10", rownames(pts_mds))

	subTitle <- paste(subTitles[expId], " ", matrixNames[expId], " (stress ", round(stress_mds, 2),")", sep = "")
	
	pts_mds$plots <- factor(pts_mds$plots,levels = unique(pts_mds$plots))
	# Plot MDS ordination
	mdsp <- ggplot(pts_mds, aes(x = MDS1, y = MDS2, color = plots)) + 
			geom_point(size = 3) + 
			geom_text(aes(label = rownames(pts_mds)), size = 3, vjust = 2) +
			geom_polygon(aes(group = pts_mds$plots), fill = NA) +
			theme(legend.position="top", legend.title=element_blank()) +
			ggtitle(subTitle) 

   assign(paste('mdsp', expId, sep=''), mdsp)
}

mdsp1 <- mdsp1 + theme(axis.title.x=element_blank(), legend.position="none")
mdsp2 <- mdsp2 + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none")
mdsp3 <- mdsp3 + theme(axis.title.x=element_blank(), legend.position="none")
mdsp4 <- mdsp4 + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none")
mdsp5 <- mdsp5 + theme(legend.position="none")
mdsp6 <- mdsp6 + theme(axis.title.y=element_blank(), legend.position="none")

pdf(paste(workingPath, "figures/mds-subplots-Bray-Curtis-", otuThr, ".pdf", sep = ""), width=7, height=10)	   
mylegend<-g_legend(mdsp)
grid.arrange(mylegend,arrangeGrob(mdsp1,mdsp2,mdsp3,mdsp4,mdsp5,mdsp6,ncol = 2, nrow=3), ncol=1, nrow=2, heights=c(1/20,19/20))
invisible(dev.off()) 
             
######## Horn-Morisita #######

for (expId in 1:n) {	
    # "-by-plot" trigger merge 2 subplots columns
    communityMatrix <- init(expId, otuThr, "-by-subplot")    
    matrixName <- matrixNames[expId]
    d.hornMorisita <- vegdist(communityMatrix, method="horn", binary=FALSE)
		
        # Run metaMDS, get points and stress
	mds <- metaMDS(d.hornMorisita)
	pts_mds <- as.data.frame(mds$points)
	stress_mds <- mds$stress
	
	# Get sample labels/factors
	pts_mds$plots <- sapply(strsplit(rownames(pts_mds), "-"), "[", 1)
	# convert plot names
	pts_mds$plots <- gsub("CM30C30", "Plot9", pts_mds$plots)
	pts_mds$plots <- gsub("LB1", "Plot10", pts_mds$plots)
	rownames(pts_mds) <- gsub("CM30C30", "Plot9", rownames(pts_mds))
	rownames(pts_mds) <- gsub("LB1", "Plot10", rownames(pts_mds))

	subTitle <- paste(subTitles[expId], matrixNames[expId], "(dissimilarity ", round(stress_mds, 2),")")
	
	pts_mds$plots <- factor(pts_mds$plots,levels = unique(pts_mds$plots))
	# Plot MDS ordination
	mdsp <- ggplot(pts_mds, aes(x = MDS1, y = MDS2, color = plots)) + 
			geom_point(size = 3) + 
			geom_text(aes(label = rownames(pts_mds)), size = 3, vjust = 2) +
			geom_polygon(aes(group = pts_mds$plots), fill = NA) +
			theme(legend.position="top", legend.title=element_blank()) +
			ggtitle(subTitle) 

   assign(paste('mdsp', expId, sep=''), mdsp)
}

mdsp1 <- mdsp1 + theme(axis.title.x=element_blank(), legend.position="none")
mdsp2 <- mdsp2 + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none")
mdsp3 <- mdsp3 + theme(axis.title.x=element_blank(), legend.position="none")
mdsp4 <- mdsp4 + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none")
mdsp5 <- mdsp5 + theme(legend.position="none")
mdsp6 <- mdsp6 + theme(axis.title.y=element_blank(), legend.position="none")

pdf(paste(workingPath, "figures/mds-subplots-Horn-Morisita-", otuThr, ".pdf", sep = ""), width=6, height=9)	
mylegend<-g_legend(mdsp)
grid.arrange(mylegend,arrangeGrob(mdsp1,mdsp2,mdsp3,mdsp4,mdsp5,mdsp6,ncol = 2, nrow=3), ncol=1, nrow=2, heights=c(1/20,19/20))
invisible(dev.off()) 
