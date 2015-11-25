
if(!exists("matrixNames")) stop("matrix names are missing !")

# source("Modules/init.r")
source("Modules/phylorare.R")

# min.size: skip if rare.min < min.size
createPhylorareTable <- function(expId, isPlot, rmSingleton, taxa.group, treeFileStem, min.size=100, verbose=T) {
  phyloTree <- getPhyloTree(treeFileStem)
  
  if (is.null(phyloTree)) {
    cat("\nWarning: no tree file, skip", taxa.group, "subset from", matrixNames[expId], ".\n")
    return(FALSE)
  }
  
  communityMatrix <- getCommunityMatrixT(expId, isPlot, rmSingleton, taxa.group)
  
  rare.max <- max(rowSums(communityMatrix))
  rare.min <- min(rowSums(communityMatrix))
  
  if (rare.min < min.size) {
    cat("\nWarning: min sample size", rare.min, "<", min.size, ", skip", taxa.group, "subset from", matrixNames[expId], ".\n")
    return(FALSE)
  }
  sample.sizes <- c(round(exp(seq(log(1), log(rare.min), length.out = 5)), digits = 0), 
                    round(exp(seq(log(rare.min), log(rare.max), length.out = 10)), digits = 0))
  
  for (ss in sample.sizes) {
    if (verbose)
      cat("Phylo rare:", matrixNames[expId], "sample size =", ss, ", subsampling by individual.\n") 
    
    # individual (default), site or species
    phylo.rare <- phylorare(communityMatrix, phyloTree, m=ss, subsampling = "individual", replace =F)
    
    if (which(sample.sizes == ss) == 1) 
      phylo.rare.df <- data.frame(row.names=rownames(phylo.rare), check.names=FALSE)
    
    if (! all(tolower(rownames(phylo.rare.df)) == tolower(rownames(phylo.rare))) )
      stop("Sample names do not match between phylo.rare and community matrix !")
    
    phylo.rare.df[,paste("size.", ss, sep="")] <- phylo.rare[,1]
  }
  
  outputFile <- file.path(workingPath, "data", paste(matrixNames[expId], 
            postfix(taxa.group, isPlot, rmSingleton, sep="-"), "phylorare", "table.csv", sep = "-"))
  write.csv(phylo.rare.df, outputFile, row.names=TRUE, quote=FALSE)
  
  return(TRUE)
}

# pathFileStem = filePath + fileStem
plotPhylorareTable <- function(phylo.rare, env, pathFileStem, verbose=T) {
  merge.df <- merge(phylo.rare.df, env[,c("Elevation","ForestType")], by = "row.names", sort = FALSE)
  # rename Row.names
  colnames(merge.df)[1] <- "Samples"
  
  melt.df <- melt(merge.df, id=c("Samples","Elevation","ForestType"))
  melt.df <- melt.df[complete.cases(melt.df),]
  # remove prefix size.
  melt.df$variable  <- gsub("^.*?\\.","",melt.df$variable) 
  
  melt.df$variable <- as.numeric(melt.df$variable)
  melt.df$Samples = factor(melt.df$Samples, unique(melt.df$Samples))
  #sort(unique(melt.df$ForestType))
  melt.df$ForestType = factor(melt.df$ForestType, c("VS2","VS3","VS5","WF7","WF9","WF11","WF12","WF13","MF20","unknown"))
  
  end.lines <- aggregate(cbind(variable, value) ~ Samples, melt.df, max)
  end.lines <- merge(end.lines, unique(melt.df[,c("Samples","Elevation","ForestType")]), by = "Samples", sort = FALSE)
  end.lines$ForestType = factor(end.lines$ForestType, c("VS2","VS3","VS5","WF7","WF9","WF11","WF12","WF13","MF20","unknown"))
  
  p <- ggplot(data=melt.df, aes(x=variable, y=value, group=Samples, colour=Elevation)) + 
    scale_shape_manual(values=c(15,16,17,0,1,2,5,6,3,4)) +
    scale_colour_gradientn(colours = c("blue", "orange")) +
    scale_x_sqrt(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) + 
    #scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), label=scientific_10) +
    xlab("Rarefaction Size") + ylab("Mean of Rooted PD") + 
    ggtitle(paste(matrixNames[expId], taxa.group)) +     
    geom_line(size=.5, alpha=0.75) +     
    #geom_point(data=end.lines, aes(shape=ForestType, colour=Elevation), size=3) +
    geom_text(data=end.lines, aes(label=Samples), size=3, alpha=0.5, hjust=-0.1, vjust=-0.6) + #, position=position_jitter(width=-1, height=3)) + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank()) 
  
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  
  pdf(paste(pathFileStem, "pdf", sep = "."), width=8, height=6)
  #print(p)
  print(grid.draw(gt))
  invisible(dev.off())   
}
