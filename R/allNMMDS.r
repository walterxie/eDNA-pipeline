
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
legend_nrow = 1

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
  env[,shapedBy] <-  substr(env[,shapedBy], 1, 2)
  n_shape <- length(unique(env[,shapedBy]))
  myshape <- seq(0, (0 + n_shape-1))
  
  env
}

######## MDS #######

for (expId in 1:n) {	
  # isPlot, rmSingleton, taxa.group, are fixed in init, when expId == n
  diss <- getDissimilarityMatrix(expId, isPlot, rmSingleton, diss.fun, taxa.group)
  
  # Run metaMDS, get points and stress
  mds <- metaMDS(dist(diss))
  pts_mds <- as.data.frame(mds$points)
  pts_mds <- pts_mds[order(rownames(pts_mds)),]
  stress_mds <- mds$stress
  
  if (expId == n) {
    linkedBy = NULL
    env <- getEnvBy(TRUE, linkedBy, colouredBy, shapedBy)
    fname <- paste("nmmds", matrixNames[expId], postfix("all", TRUE, FALSE, sep="-"), diss.fun, sep = "-")
  } else {
    env <- getEnvBy(FALSE, linkedBy, colouredBy, shapedBy)
    fname <- paste("nmmds", matrixNames[expId], postfix(taxa.group, TRUE, rmSingleton, sep="-"), diss.fun, sep = "-")
  }
  
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

  # Plot MDS ordination
  mdsp <- ggplot(pts_mds_env, aes_string(x="MDS1", y="MDS2", color=colouredBy)) + 
    geom_text(aes(label=Row.names), size = 3, vjust = 2) +
    geom_point(aes_string(shape=shapedBy), size = 3) + scale_shape(solid=FALSE) +
    theme_bw() +
    theme(legend.position="top", panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank()) +
    guides(col=guide_legend(nrow=legend_nrow)) +
    ggtitle(plotTitle) 

  if (!is.null(linkedBy)) 		
    mdsp <- mdsp + geom_polygon(data = hulls, aes_string(group=linkedBy), fill = NA) 		
  
  gt <- ggplot_gtable(ggplot_build(mdsp))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  
  pdf(paste(workingPath, figDir, "/", fname, ".pdf", sep = ""), width=5, height=6, useDingbats = FALSE) 
  print(grid.draw(gt))
  invisible(dev.off()) 
}
