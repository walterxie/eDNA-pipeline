
library(ggplot2)
library(vegan)
library(vegetarian)
library(grid)
library(gridExtra)

if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("isPlot")) isPlot = FALSE # by subplot
if(!exists("otuThr")) otuThr = 97
if(!exists("diss.fun")) diss.fun="beta1-1"
if(!exists("taxa.group")) taxa.group="all"


n <- length(matrixNames) 
subTitles <- c("a) ","b) ","c) ","d) ","e) ","f) ")

source("Modules/init.R", local=TRUE)

cat("Analysis: Non-metric MDS of", taxa.group, "taxa group(s) using", diss.fun, "dissimilarity. \n")

######## env #######
linkedBy = "Plot"
colouredBy = "Elevation"
shapedBy = "ForestType"  

cat("isPlot =", isPlot, "linkedBy =", linkedBy, ", colouredBy =", colouredBy, ", shapedBy =", shapedBy, ".\n")

getEnvBy <- function(isPlot, linkedBy, colouredBy, shapedBy) {
  env <- getSampleMetaData(isPlot)
  
  if (! colouredBy %in% colnames(env))
    stop("Error: colouredBy is not defined from colnames of env  !")
  if (! is.null(linkedBy) && ! linkedBy %in% colnames(env) )
    stop("Error: linkedBy is not defined from colnames of env  !")
  if (! shapedBy %in% colnames(env) )
    stop("Error: shapedBy is not defined from colnames of env  !")
  
  env[,shapedBy] <- gsub(":.*", "", env[,shapedBy], ignore.case = T)
  #env[,shapedBy] <-  substr(env[,shapedBy], 1, 2)
  
  env
}

######## MDS #######

for (expId in 1:n) {	
  # isPlot, rmSingleton, taxa.group, are fixed in init, when expId == n
  diss <- getDissimilarityMatrix(expId, FALSE, rmSingleton, diss.fun, taxa.group)
  
  # Run metaMDS, get points and stress
  mds <- metaMDS(as.dist(diss))
  pts_mds <- as.data.frame(mds$points)
  pts_mds <- pts_mds[order(rownames(pts_mds)),]
  stress_mds <- mds$stress
  
  if (expId == n) {
    linkedBy = NULL
    env <- getEnvBy(TRUE, linkedBy, colouredBy, shapedBy)
    fname <- paste("nmmds", matrixNames[expId], postfix("all", TRUE, FALSE, sep="-"), diss.fun, sep = "-")
  } else {
    env <- getEnvBy(FALSE, linkedBy, colouredBy, shapedBy)
    fname <- paste("nmmds", matrixNames[expId], postfix(taxa.group, FALSE, rmSingleton, sep="-"), diss.fun, sep = "-")
  }
  rownames(pts_mds) <- toupper(rownames(pts_mds))
  rownames(env) <- toupper(rownames(env))
  
  if (is.null(linkedBy)) {
    pts_mds_env <- merge(pts_mds, env[,c(colouredBy, shapedBy)], by = "row.names")
  } else {
    pts_mds_env <- merge(pts_mds, env[,c(linkedBy, colouredBy, shapedBy)], by = "row.names")
  }
  
  if (nrow(pts_mds_env) != nrow(pts_mds)) 
    stop("Error: sample names in community matrix must be a subset of sample names in environmental data file !")
  
  if (!is.null(linkedBy)) {
    # Convex hull http://stackoverflow.com/questions/16428962/convex-hull-ggplot-using-data-tables-in-r
    pts_mds_dt <- data.table(pts_mds_env, key = linkedBy)
    hulls <- pts_mds_dt[, .SD[chull(MDS1, MDS2)], by = linkedBy]
  }

  # subTitles[expId], 
  plotTitle <- paste(matrixNames[expId], " (stress ", round(stress_mds, 2), ")", sep = "")
  
  n_shape <- length(unique(pts_mds_env[,shapedBy]))
  myshape <- seq(5, (5 + n_shape-1))
  
  # Plot MDS ordination
  mdsp <- ggplot(pts_mds_env, aes_string(x="MDS1", y="MDS2")) + 
    geom_point(aes_string(shape=shapedBy, color=colouredBy), size = 3) + 
    scale_shape_manual(values=myshape) + # The shape palette can deal with a maximum of 6 discrete values
    geom_text(aes_string(label="Row.names", color=colouredBy), size = 3, vjust = 2, alpha = 0.5) +
    theme_bw() + scale_colour_gradientn(colours = c("red", "green", "blue")) +
    theme(legend.position="top", plot.title = element_text(size = 8), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank()) +
    guides(shape=guide_legend(nrow=2)) +
    ggtitle(plotTitle) 

  if (!is.null(linkedBy)) 		
    mdsp <- mdsp + geom_polygon(data = hulls, aes_string(mapping=linkedBy, color=colouredBy), fill = NA) 		
  
  gt <- ggplot_gtable(ggplot_build(mdsp))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  
  pdf(paste(workingPath, figDir, "/", fname, ".pdf", sep = ""), width=5, height=6, useDingbats = FALSE) 
  print(grid.draw(gt))
  invisible(dev.off()) 
}
