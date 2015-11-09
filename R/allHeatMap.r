# http://stackoverflow.com/questions/6673162/reproducing-lattice-dendrogram-graph-with-ggplot2
# improved by Walter Xie
# 2 Nov 2015

library(ggplot2)
library(reshape2)
library(ggdendro)
library(plyr)

if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("diss.fun")) diss.fun="beta1-1"
if(!exists("taxa.group")) taxa.group="assigned"

n <- length(matrixNames) 

source("Modules/init.R", local=TRUE)

### Set up a blank theme
theme_none <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(colour=NA),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(), 
  axis.line = element_blank(),
  axis.ticks.length = unit(0,"null")
)

# main
n_cluster=9
mypalette <- rainbow(n_cluster)
cat("Analysis: Heat-map with dendrogram on", n_cluster, "clusters. \n")

for (expId in 1:n) {
  ### eDNA ###
  isP <- isPlot
  min2 <- rmSingleton
  most.abundant=100
  x.lable=paste(most.abundant, "most abundant OTUs")
  bar.label="Reads"
  
  if (expId == n) {
    ### Veg ###
    isP <- TRUE
    min2 <- FALSE
    x.lable="Species"
    bar.label="SBA"
  } 
  
  communityMatrix <- getCommunityMatrixT(expId, isP, min2, taxa.group)
  # arrange the species from most abundant on the left to least abundant on the right
  communityMatrix <- communityMatrix[, order(colSums(communityMatrix), decreasing = TRUE)]
  
  if (expId < n) {
    communityMatrix <- communityMatrix[, 1:most.abundant]
    
    taxaPaths <- getTaxaPaths(expId, taxa.group)
    
    taxaAssg <- merge(t(communityMatrix), taxaPaths, by = "row.names", sort = FALSE)
    # move 1st col Row.names to row.names
    rownames(taxaAssg) <- taxaAssg[,"Row.names"]
    taxaAssg <- taxaAssg[,-1]
    
    if (! all(colnames(communityMatrix) == rownames(taxaAssg)) )
      stop("OTUs names in community matrix do not match taxa classifications !")
    
    taxaAssg[,"phylum"] <- gsub("(\\s\\[=.*\\])", "", taxaAssg[,"phylum"])
    taxaAssg[,"class"] <- gsub("(\\s\\[=.*\\])", "", taxaAssg[,"class"])
    
    colnames(communityMatrix) <- paste(taxaAssg[,"kingdom"], taxaAssg[,"phylum"], taxaAssg[,"class"], 1:nrow(taxaAssg), sep="_")
    
    #colnames(communityMatrix) <- paste(sapply(strsplit(colnames(communityMatrix), "_"), "[[", 6), 
    #                                   sapply(strsplit(colnames(communityMatrix), "_"), "[[", 7), sep="_")
    #colnames(communityMatrix) <- gsub("\\|gene", "", colnames(communityMatrix), ignore.case = T)
    #colnames(communityMatrix) <- paste("OTU", colnames(communityMatrix), sep="_")
    
    #if (length(unique(colnames(communityMatrix))) != ncol(communityMatrix)) 
    #  stop("Find duplication from simplified OTU names")
  }
  
  max.cm <- max(communityMatrix)
  my_breaks <- round_any(exp(seq(log(2), log(max.cm), length = 5)), 10)
  my_breaks[my_breaks==0] <- 1
  
  ### Dendrogram plot
  # isPlot, rmSingleton, taxa.group, are fixed in init, when expId == n
  diss <- getDissimilarityMatrix(expId, isP, min2, diss.fun, taxa.group)
  
  hc <- hclust(as.dist(diss), "ave") # "average" (= UPGMA)
  
  # cut dendrogram in n_cluster clusters
  clusMember = cutree(hc, n_cluster)
  # using dendrogram objects
  hcd = as.dendrogram(hc)
  
  dd.col <- as.dendrogram(hcd)
  col.ord <- order.dendrogram(dd.col)
  
  ddata_y <- dendro_data(dd.col)
  
  # arrange the plots the same way as dendrograms labels
  communityMatrix <- communityMatrix[as.vector(ddata_y$labels$label),]
  
  ### Create data for plot ###
  
  communityMatrix$Samples = rownames(communityMatrix)
  mdf <- melt(communityMatrix)
  mdf$Samples <- factor(mdf$Samples, levels=unique(mdf$Samples))
  
  ### Create plot components ###    
  # na.value = "#132B43", high = "#56B1F7", low = "#132B43",
  # Heatmap
  p1 <- ggplot(mdf, aes(x=variable, y=Samples)) + 
    geom_tile(aes(fill=value)) + 
    scale_fill_gradient(na.value = "transparent", low = "#56B1F7", high = "#132B43",
                        name = bar.label, trans = "log", breaks = my_breaks, labels = my_breaks) +
    #guides(fill = guide_legend(reverse=TRUE)) +
    xlab(x.lable) +
    theme(axis.title.y=element_blank(), axis.text.x=element_text(angle=60, hjust=1), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
  
  # naive method to flip over tree
  ddata_y$segments$yend <- max(ddata_y$segments$y) + min(ddata_y$segments$y) - ddata_y$segments$yend
  ddata_y$segments$y <- max(ddata_y$segments$y) + min(ddata_y$segments$y) - ddata_y$segments$y
  
  # Dendrogram plot
  p2 <- ggplot(segment(ddata_y)) + 
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
    coord_flip() + theme_none
  
  ### Draw graphic ###
  if (expId == n) {
    fname <- paste("heatmap", matrixNames[expId], postfix("all", isP, min2, sep="-"), diss.fun, sep = "-")
    pdf.width=15 
    pdf.height=6
    p2.vp=viewport(0.1, 0.95, x=0.06, y=0.55)
    p1.vp=viewport(0.92, 1, x=0.55, y=0.5)
  } else if (isP) {
    fname <- paste("heatmap", matrixNames[expId], postfix(taxa.group, isP, min2, sep="-"), diss.fun, sep = "-")
    pdf.width=18 
    pdf.height=7
    p2.vp=viewport(0.1, 0.91, x=0.06, y=0.57)
    p1.vp=viewport(0.92, 1, x=0.55, y=0.5)
  } else {
    fname <- paste("heatmap", matrixNames[expId], postfix(taxa.group, isP, min2, sep="-"), diss.fun, sep = "-")
    pdf.width=18 
    pdf.height=10
    p2.vp=viewport(0.1, 0.8, x=0.06, y=0.63)
    p1.vp=viewport(0.92, 1, x=0.55, y=0.5)
  }
   
  pdf(paste(workingPath, figDir, "/", fname, ".pdf", sep = ""), width=pdf.width, height=pdf.height)
  
  grid.newpage()
  print(p2, vp=p2.vp)
  print(p1, vp=p1.vp)
  
  invisible(dev.off())
}
