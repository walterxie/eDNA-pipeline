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

cat("Analysis: eDNA community phylogenetic structure. \n")

for (expId in 1:(n-1)) {
    ### eDNA ###
    communityMatrix <- getCommunityMatrixT(expId, isPlot, rmSingleton, taxa.group)
    
    phyloTree <- getPhyloTree(treeFile)

    cm.env <- getSampleMetaData(isPlot)  
    
    comm.phylo.struc(communityMatrix, phyloTree, treeFileStem, tableFile, verbose)
}
