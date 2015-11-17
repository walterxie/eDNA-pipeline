# phylorare: http://davidnipperess.blogspot.co.nz/search/label/Software


if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("tableFile")) stop("table file is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("taxa.group")) taxa.group="assigned"
if(!exists("isPlot")) isPlot = FALSE # by subplot


n <- length(matrixNames) 

source("Modules/init.r")
source("Modules/phylorare.R", local=TRUE)

cat("Intermediate data: create rarefaction rarefied Phylogenetic Diversity table",
    ", rmSingleton =", rmSingleton, ", isPlot =", isPlot, ", otuThr =", otuThr, "\n") 

for (expId in 1:(n-1)) {
    ### eDNA ###
    # 16S-assigned-min2
    treeFileStem <- tolower(paste(matrixNames[expId], postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-"))
    phyloTree <- getPhyloTree(treeFileStem)

    if (is.null(phyloTree)) {
      cat("\nSkip", taxa.group, "subset from", matrixNames[expId], "because of no tree file.\n")
      next
    }
    
    communityMatrix <- getCommunityMatrixT(expId, isPlot, rmSingleton, taxa.group)
    
    rare.max <- max(rowSums(communityMatrix))
    sample.sizes <- round(seq(1, rare.max, length.out = 11), digits = 0)
    
    for (ss in sample.sizes) {
      if (verbose)
        cat("Phylo rare:", matrixNames[expId], "sample size =", ss, "\n") 
      
      # individual (default), site or species
      phylo.rare <- phylorare(communityMatrix, phyloTree, m=ss, subsampling = "individual", replace =F)
      
      if (ss == 1) 
        phylo.rare.df <- data.frame(row.names=rownames(phylo.rare), check.names=FALSE)
       
      if (! all(tolower(rownames(phylo.rare.df)) == tolower(rownames(phylo.rare))) )
        stop("Sample names do not match between phylo.rare and community matrix !")
      
      phylo.rare.df[,paste("size.", ss, sep="")] <- phylo.rare[,1]
    }
    
    outputRFTable <- file.path(workingPath, "data", paste(matrixNames[expId], 
          postfix(taxa.group, isPlot, rmSingleton, sep="-"), "phylorare", "table.csv", sep = "-"))
    write.csv(phylo.rare.df, outputRFTable, row.names=TRUE, quote=FALSE)
}
