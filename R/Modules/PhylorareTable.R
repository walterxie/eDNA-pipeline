
if(!exists("matrixNames")) stop("matrix names are missing !")

# source("Modules/init.r")
source("Modules/phylorare.R")

# min.size: skip if rare.min < min.size
createPhylorareTable <- function(expId, isPlot, rmSingleton, taxa.group, treeFileStem, min.size=100, verbose=T) {
  pdFilePath <- file.path(workingPath, "data", "pdrf")
  mkdir(pdFilePath) 
  
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
  sample.sizes <- c(round(exp(seq(log(1), log(rare.min), length.out = 6)), digits = 0)[1:5], 
                    round(exp(seq(log(rare.min), log(rare.max), length.out = 9)), digits = 0))
  
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
  
  outputFile <- file.path(pdFilePath, paste(matrixNames[expId], 
            postfix(taxa.group, isPlot, rmSingleton, sep="-"), "phylorare", "table.csv", sep = "-"))
  write.csv(phylo.rare.df, outputFile, row.names=TRUE, quote=FALSE)
  
  return(TRUE)
}

