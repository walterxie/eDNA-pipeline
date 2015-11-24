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
source("Modules/PhylorareTable.R")

cat("Intermediate data: create rarefaction rarefied Phylogenetic Diversity table",
    ", rmSingleton =", rmSingleton, ", isPlot =", isPlot, ", otuThr =", otuThr, "\n") 

for (expId in 1:(n-1)) {
    ### eDNA ###
    # 16S-assigned-min2, isPlot always FALSE for a valid treeFileStem
    treeFileStem <- tolower(paste(matrixNames[expId], postfix(taxa.group, FALSE, rmSingleton, sep="-"), sep = "-"))
    
    createPhylorareTable(expId, isPlot, rmSingleton, taxa.group, treeFileStem, min.size=100, verbose=verbose)
}
