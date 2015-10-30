# http://stackoverflow.com/questions/6673162/reproducing-lattice-dendrogram-graph-with-ggplot2

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
  axis.text.y = element_blank(), #element_text(angle=90, hjust=1),
  axis.line = element_blank()
  #axis.ticks.length = element_blank()
)

# main
for (expId in 1:n) {
  ### eDNA ###
  isP <- isPlot
  min2 <- rmSingleton
  inds="OTUs"
  
  if (expId == n) {
    ### Veg ###
    isP <- TRUE
    min2 <- FALSE
    taxa.group="all"
    inds="Species"
  } 
  
  communityMatrix <- getCommunityMatrixT(expId, isP, min2, taxa.group)
  max.cm <- max(communityMatrix)
  my_breaks <- round_any(exp(seq(log(2), log(max.cm), length = 5)), 10)
  my_breaks[my_breaks==0] <- 1
  
  ### Dendrogram plot
  # isPlot, rmSingleton, taxa.group, are fixed in init, when expId == n
  diss <- getDissimilarityMatrix(expId, isP, min2, diss.fun, taxa.group)
  
  hc <- hclust(as.dist(diss), "ave") # "average" (= UPGMA)
  
  dd.col <- as.dendrogram(hc)
  col.ord <- order.dendrogram(dd.col)
  
  ddata_y <- dendro_data(dd.col)
  
  ### Create data for plot ###
  
  communityMatrix$Samples = rownames(communityMatrix)
  mdf <- melt(communityMatrix)
  
  ### Create plot components ###    
  # Heatmap
  p1 <- ggplot(mdf, aes(x=variable, y=Samples)) + 
    geom_tile(aes(fill=value)) + 
    scale_fill_gradient(na.value = "transparent", low = "#56B1F7", high = "#132B43",
                        name = "SBA", trans = "log", breaks = my_breaks, labels = my_breaks) +
    #guides(fill = guide_legend(reverse=TRUE)) +
    theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
  
  # Dendrogram plot
  p2 <- ggplot(segment(ddata_y)) + 
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
    coord_flip() + theme_none
  
  ### Draw graphic ###
  fname <- paste("heatmap", matrixNames[expId], postfix(taxa.group, isP, min2, sep="-"), diss.fun, sep = "-")
  pdf(paste(workingPath, figDir, "/", fname, ".pdf", sep = ""), width=12, height=6)
  
  grid.newpage()
  print(p1, vp=viewport(0.91, 1, x=0.45, y=0.5))
  print(p2, vp=viewport(0.09, 1, x=0.93, y=0.55))
  
  invisible(dev.off())
}
