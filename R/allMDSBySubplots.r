
library(ggplot2)
library(vegan)
library(vegetarian)
library(grid)
library(gridExtra)

# change config below
#figDir <- "figures"
#sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
#setwd(sourcePath)
#workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
#matrixNames <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name   

if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")
if(!exists("levels")) stop("levels of Jost diversity are missing !")
if(!exists("qs")) stop("qs of Jost diversity are missing !")

subTitles <- c("(a)","(b)","(c)","(d)","(e)","(f)")

n <- length(matrixNames) 
mypalette <- c("red", "orange", "green", "purple", "blue", "brown")
myshape <- seq(0, (0 + n-1))

otuThr = 97

source("Modules/init.R", local=TRUE)
  
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
	pts_mds$plots <- gsub("CM30C30", "Plot9", pts_mds$plots, ignore.case = T)
	pts_mds$plots <- gsub("LB1", "Plot10", pts_mds$plots, ignore.case = T)
	rownames(pts_mds) <- gsub("CM30C30", "Plot9", rownames(pts_mds), ignore.case = T)
	rownames(pts_mds) <- gsub("LB1", "Plot10", rownames(pts_mds), ignore.case = T)
	rownames(pts_mds) <- gsub("[Plot]", "", rownames(pts_mds))

	subTitle <- paste(subTitles[expId], " ", matrixNames[expId], " (stress ", round(stress_mds, 2),")", sep = "")
	
	pts_mds$plots <- factor(pts_mds$plots,levels = unique(pts_mds$plots))
	# Plot MDS ordination
	# both position_jitter and directlabels not working
	mdsp <- ggplot(pts_mds, aes(x = MDS1, y = MDS2, color = plots)) + 
			geom_point(size = 3) + 
			geom_text(aes(label = rownames(pts_mds)), size = 3, vjust = 2) +
			geom_polygon(aes(group = pts_mds$plots), fill = NA) + theme_bw() +
			theme(legend.position="top", legend.title=element_blank(), panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(), panel.background = element_blank()) +
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

for (expId in 1:n) {
	mdsp <- get(paste('mdsp', expId, sep=''))
	gg2 <- ggplot_gtable(ggplot_build(mdsp))
	gg2$layout$clip[gg2$layout$name == "panel"] <- "off"
	#grid.draw(gg2)
	assign(paste('mdsp', expId, sep=''), gg2)
}

pdf(paste(workingPath, figDir, "/mds-subplots-beta1-", otuThr, ".pdf", sep = ""), width=7, height=10, useDingbats = FALSE)	   
grid.arrange(mylegend,arrangeGrob(mdsp1,mdsp2,mdsp3,mdsp4,mdsp5,mdsp6,ncol = 2, nrow=3), ncol=1, nrow=2, heights=c(1/20,19/20))
invisible(dev.off()) 

print(paste("Plot mds-subplots-beta1-", otuThr, ".pdf to folder : ", workingPath, figDir, sep=""))  

######## Jaccard #######

for (expId in 1:n) {	
    # "-by-plot" trigger merge 2 subplots columns
    communityMatrix <- init(expId, otuThr, "-by-subplot")    
    matrixName <- matrixNames[expId]
    d.jaccard <- vegdist(communityMatrix, method="jaccard")
		
    # Run metaMDS, get points and stress
	mds <- metaMDS(d.jaccard)
	pts_mds <- as.data.frame(mds$points)
	stress_mds <- mds$stress
	
	# Get sample labels/factors
	pts_mds$plots <- sapply(strsplit(rownames(pts_mds), "-"), "[", 1)
	# convert plot names
	pts_mds$plots <- gsub("CM30C30", "Plot9", pts_mds$plots, ignore.case = T)
	pts_mds$plots <- gsub("LB1", "Plot10", pts_mds$plots, ignore.case = T)
	rownames(pts_mds) <- gsub("CM30C30", "Plot9", rownames(pts_mds), ignore.case = T)
	rownames(pts_mds) <- gsub("LB1", "Plot10", rownames(pts_mds), ignore.case = T)
	rownames(pts_mds) <- gsub("[Plot]", "", rownames(pts_mds))

	subTitle <- paste(subTitles[expId], " ", matrixNames[expId], " (stress ", round(stress_mds, 2),")", sep = "")
	
	pts_mds$plots <- factor(pts_mds$plots,levels = unique(pts_mds$plots))
	# Plot MDS ordination
	mdsp <- ggplot(pts_mds, aes(x = MDS1, y = MDS2, color = plots)) + 
			geom_point(size = 3) + 
			geom_text(aes(label = rownames(pts_mds)), size = 3, vjust = 2) +
			geom_polygon(aes(group = pts_mds$plots), fill = NA) + theme_bw() +
			theme(legend.position="top", legend.title=element_blank(), panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(), panel.background = element_blank()) +
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

for (expId in 1:n) {
	mdsp <- get(paste('mdsp', expId, sep=''))
	gg2 <- ggplot_gtable(ggplot_build(mdsp))
	gg2$layout$clip[gg2$layout$name == "panel"] <- "off"
	#grid.draw(gg2)
	assign(paste('mdsp', expId, sep=''), gg2)
}

pdf(paste(workingPath, figDir, "/mds-subplots-jaccard-", otuThr, ".pdf", sep = ""), width=7, height=10, useDingbats = FALSE)	   
grid.arrange(mylegend,arrangeGrob(mdsp1,mdsp2,mdsp3,mdsp4,mdsp5,mdsp6,ncol = 2, nrow=3), ncol=1, nrow=2, heights=c(1/20,19/20))
invisible(dev.off()) 
             
print(paste("Plot mds-subplots-jaccard-", otuThr, ".pdf to folder : ", workingPath, figDir, sep="")) 

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
	pts_mds$plots <- gsub("CM30C30", "Plot9", pts_mds$plots, ignore.case = T)
	pts_mds$plots <- gsub("LB1", "Plot10", pts_mds$plots, ignore.case = T)
	rownames(pts_mds) <- gsub("CM30C30", "Plot9", rownames(pts_mds), ignore.case = T)
	rownames(pts_mds) <- gsub("LB1", "Plot10", rownames(pts_mds), ignore.case = T)
	rownames(pts_mds) <- gsub("[Plot]", "", rownames(pts_mds))

	subTitle <- paste(subTitles[expId], matrixNames[expId], "(stress ", round(stress_mds, 2),")")
	
	pts_mds$plots <- factor(pts_mds$plots,levels = unique(pts_mds$plots))
	# Plot MDS ordination
	mdsp <- ggplot(pts_mds, aes(x = MDS1, y = MDS2, color = plots)) + 
			geom_point(size = 3) + 
			geom_text(aes(label = rownames(pts_mds)), size = 3, vjust = 2) +
			geom_polygon(aes(group = pts_mds$plots), fill = NA) + theme_bw() +
			theme(legend.position="top", legend.title=element_blank(), panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(), panel.background = element_blank()) +
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

for (expId in 1:n) {
	mdsp <- get(paste('mdsp', expId, sep=''))
	gg2 <- ggplot_gtable(ggplot_build(mdsp))
	gg2$layout$clip[gg2$layout$name == "panel"] <- "off"
	#grid.draw(gg2)
	assign(paste('mdsp', expId, sep=''), gg2)
}

pdf(paste(workingPath, figDir, "/mds-subplots-Horn-Morisita-", otuThr, ".pdf", sep = ""), width=7, height=10, useDingbats = FALSE)	
grid.arrange(mylegend,arrangeGrob(mdsp1,mdsp2,mdsp3,mdsp4,mdsp5,mdsp6,ncol = 2, nrow=3), ncol=1, nrow=2, heights=c(1/20,19/20))
invisible(dev.off()) 

print(paste("Plot mds-subplots-Horn-Morisita-", otuThr, ".pdf to folder : ", workingPath, figDir, sep="")) 
