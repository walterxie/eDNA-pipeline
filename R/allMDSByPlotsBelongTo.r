



library(ggplot2)
library(vegan)
library(vegetarian)
library(grid)
library(gridExtra)
library(data.table)

####### turn over = beta 1 -1 ########
turnover <- function(communityMatrix) {    
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
#subTitles <- c("(a)","(b)","(c)","(d)","(e)","(f)")
taxaFiles <- c("DOCx9_16S_R1_reads_OTUs_num_nt_taxatable.txt", "DOCx9_18S_PEARcontigs_ME1.0_OTUs_num_nt_taxatable.txt", 
				"DOCx9_26S_R1_reads_OTUs_num_nt_taxatable.txt", "DOCx9_ITS_PEARcontigs_ME1.0_OTUs_num_nt_taxatable.txt", 
				"DOCx9_FolCO1_R1_reads_OTUs_num_nt_taxatable.txt", "DOCx9_ShCO1_PEARcontigs_ME1.0_OTUs_num_nt_taxatable.txt")

n <- length(matrixNames) 
mypalette <- c("red", "orange", "green", "purple", "blue", "brown")
myshape <- seq(0, (0 + n-1))

otuThr = 97


# not make a plot for the data set whose reads <= minReads
minReads = 100

#belongTo = "Insecta" 
belongToList = c("Bacteria", "Fungi", "Viridiplantae", "Arthropoda", "Metazoa", "Insecta") # they are in different rank

linkedBy = "Vegetation"


######## environmental data #######
elev <- read.table(paste(workingPath, "data/plot_elevations.txt", sep=""), sep="\t", header=T, row.names=1)
elev <- elev[order(rownames(elev)), ]
elev <- elev[-which(rownames(elev)=="AL132"),]


for (belongId in 1:length(belongToList)) {
	belongTo = belongToList[belongId]
	cat(paste("\nChoose tax group ", belongTo, " :\n"))

	for (expId in 1:n) {	
		cat(paste("\nLoad data for", matrixNames[expId], " :\n"))

		##### load community matrix #####
		inputCM <- paste(workingPath, "data/", matrixNames[expId], ".txt", sep="") 
		communityMatrix <- read.table(inputCM, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
		communityMatrix <- communityMatrix[order(rownames(communityMatrix)),]
		communityMatrix <- communityMatrix[,order(colnames(communityMatrix))]
		# U115 = U155 (correct name) and OJ95 = CJ95 (correct name). AL132 is the "difficult" site extracting DNA 
		colnames(communityMatrix) <- gsub("OJ95", "CJ95", colnames(communityMatrix), ignore.case = T)
		colnames(communityMatrix) <- gsub("U115", "U155", colnames(communityMatrix), ignore.case = T)

		rownames(communityMatrix) <- gsub(";size.*;", "", rownames(communityMatrix))

		if ( all( tolower(rownames(elev)) != tolower(colnames(communityMatrix)) ) ) 
			stop("Site names in community matrix and environmental-elevation data file not matched !")

		##### load taxa paths #####
		inputTaxa <- paste(workingPath, "Taxonomy_tables/", taxaFiles[expId], sep="")
		taxaPaths <- read.table(inputTaxa, header=TRUE, row.names=1, sep = "\t", stringsAsFactors=FALSE)  
		taxaPaths <- taxaPaths[order(rownames(taxaPaths)),]

		rownames(taxaPaths) <- gsub(";size.*;", "", rownames(taxaPaths))

		##### filter out rows not belong to given taxa belongTo ##### 
		taxaPaths <- taxaPaths[which(grepl(belongTo, taxaPaths[,1])),] # taxaPaths[,1] is taxa path separated by ;
	
		if (belongTo == "Metazoa") { # non Arthropoda
			taxaPaths <- taxaPaths[-which(grepl("Arthropoda", taxaPaths[,1])),]
		}
	
		# not make a plot for the data set whose reads < minReads
		if (nrow(taxaPaths) < minReads) {
			print( paste("Taxa of", matrixNames[expId], "NOT have classification for", belongTo, ", jump to next in the loop.") )	
			next
		}

		##### community matrix by a taxa group #####
		communityMatrix <- communityMatrix[match(rownames(taxaPaths),rownames(communityMatrix)),]

		if ( all( tolower(rownames(communityMatrix)) != tolower(rownames(taxaPaths)) ) ) 
			stop( paste("OTU names in", matrixNames[expId], "community matrix file and taxa path file are not matched !") )

		######## beta 1 -1 #######
		# rotate 
		communityMatrix <- as.data.frame(t(communityMatrix)) 

		d.beta1 <- turnover(communityMatrix)
		
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

		subTitle <- paste(belongTo, " ", matrixNames[expId], " (stress ", round(stress_mds, 2),")", sep = "")
		
		# Convex hull http://stackoverflow.com/questions/16428962/convex-hull-ggplot-using-data-tables-in-r
		pts_mds_dt <- data.table(pts_mds,key= linkedBy)
		hulls <- pts_mds_dt[, .SD[chull(MDS1, MDS2)], by = linkedBy]

#		pts_mds[,linkedBy] <- factor(pts_mds[,linkedBy],levels = unique(pts_mds[,linkedBy]))
		# Plot MDS ordination
		mdsp <- ggplot(pts_mds, aes(x = MDS1, y = MDS2, color = Elevation)) + 
				geom_point(size = 3) + 
				geom_text(aes(label = rownames(pts_mds)), size = 3, vjust = 2) +
				geom_polygon(data = hulls, aes_string(group=linkedBy), fill = NA) +
				theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
				ggtitle(subTitle) 
				
		gt <- ggplot_gtable(ggplot_build(mdsp))
		gt$layout$clip[gt$layout$name == "panel"] <- "off"

		pdf(paste(workingPath, "figures/mds-", belongTo, "-", matrixNames[expId], ".pdf", sep = ""), width=5, height=6)	   
		print(grid.draw(gt))
		invisible(dev.off()) 
	   
	} #END for expId

} #END for belongId

