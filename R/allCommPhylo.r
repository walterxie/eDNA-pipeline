# Picante: R tools for integrating phylogenies and ecology


if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("tableFile")) stop("table file is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("taxa.group")) taxa.group="assigned"
if(!exists("isPlot")) isPlot = FALSE # by subplot


n <- length(matrixNames) 

source("Modules/CommunityPhylogenetic.R", local=TRUE)

cat("\nAnalysis: eDNA community phylogenetic structure. \n")

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
    
    cm.env <- getSampleMetaData(isPlot)  
    
    fname <- paste("pd-beta", matrixNames[expId], postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-")
    comdistFile <- file.path(workingPath, "data", paste(fname, ".csv", sep = ""))
    
    comm.phylo.struc(communityMatrix, phyloTree, treeFileStem, tableFile, verbose, comdistFile)
}
