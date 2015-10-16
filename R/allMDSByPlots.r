# "-by-plot" trigger merge 2 subplots columns
library(ggplot2)
library(vegan)
library(vegetarian)
library(grid)
library(gridExtra)
library(data.table)

######## beta 1 -1 #######
beta1 <- function(communityMatrix) {    
	# including diagonal
    rowsNum <- nrow(communityMatrix) * (nrow(communityMatrix) + 1) / 2  		
	d.beta1 <- matrix(0,nrow=nrow(communityMatrix),ncol=nrow(communityMatrix))
	colnames(d.beta1) <- c(rownames(communityMatrix))
    rownames(d.beta1) <- c(rownames(communityMatrix))
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

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


# change config below
sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/MiSeq/DOCt1/R/Modules/"
setwd(sourcePath)

workingPath <- "~/svn/compevol/research/NZGenomicObservatory/MiSeq/DOCt1/"
experiments <-  c("16S", "18S", "26S", "ITS", "FolCO1", "ShCO1") # only for cm file name and folder name   
matrixNames <-  experiments
subTitles <- c("(a)","(b)","(c)","(d)","(e)","(f)")

n <- length(matrixNames) 
mypalette <- c("red", "orange", "green", "purple", "blue", "brown")
myshape <- seq(0, (0 + n-1))

otuThr = 97

source("init.R", local=TRUE)
  
linkedBy = "Vegetation"


######## environmental data #######
elev <- read.table(paste(workingPath, "data/plot_elevations.txt", sep=""), sep="\t", header=T, row.names=1)
elev <- elev[order(rownames(elev)), ]
elev <- elev[-which(rownames(elev)=="AL132"),]


######## beta 1 -1 #######

for (expId in 1:n) {	
    communityMatrix <- initByPlot(expId) 
    d.beta1 <- beta1(communityMatrix)
		
    # Run metaMDS, get points and stress
	mds <- metaMDS(d.beta1)
	pts_mds <- as.data.frame(mds$points)
	pts_mds <- pts_mds[order(rownames(pts_mds)),]
	stress_mds <- mds$stress

	if ( all( tolower(rownames(elev)) != tolower(rownames(pts_mds)) ) ) 
		stop("Site names in MDS plots and environmental-elevation data file not matched !")

	# Get sample labels/factors
	pts_mds$Elevation <- elev$Elevation
	pts_mds[,linkedBy] <- elev[,linkedBy]

	subTitle <- paste(subTitles[expId], " ", matrixNames[expId], " (stress ", round(stress_mds, 2),")", sep = "")
	
	# Convex hull http://stackoverflow.com/questions/16428962/convex-hull-ggplot-using-data-tables-in-r
	pts_mds_dt <- data.table(pts_mds,key= linkedBy)
	hulls <- pts_mds_dt[, .SD[chull(MDS1, MDS2)], by = linkedBy]
	
#	pts_mds[,linkedBy] <- factor(pts_mds[,linkedBy],levels = unique(pts_mds[,linkedBy]))
	# Plot MDS ordination
	mdsp <- ggplot(pts_mds, aes(x = MDS1, y = MDS2, color = Elevation)) + 
			geom_text(aes(label = rownames(pts_mds)), size = 3, vjust = 2) +
			geom_point(size = 3) + 
			geom_polygon(data = hulls, aes_string(group=linkedBy), fill = NA) +
			theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
			ggtitle(subTitle) 

   assign(paste('mdsp', expId, sep=''), mdsp)
}

mylegend<-g_legend(mdsp)

mdsp1 <- mdsp1 + theme(axis.title.x=element_blank(), legend.position="none")
mdsp2 <- mdsp2 + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none")
mdsp3 <- mdsp3 + theme(axis.title.x=element_blank(), legend.position="none")
mdsp4 <- mdsp4 + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none")
mdsp5 <- mdsp5 + theme(legend.position="none")
mdsp6 <- mdsp6 + theme(axis.title.y=element_blank(), legend.position="none")

pdf(paste(workingPath, "figures/mds-beta1-", otuThr, ".pdf", sep = ""), width=7, height=10)	   
grid.arrange(mylegend,arrangeGrob(mdsp1,mdsp2,mdsp3,mdsp4,mdsp5,mdsp6,ncol = 2, nrow=3), ncol=1, nrow=2, heights=c(1/20,19/20))
invisible(dev.off()) 


