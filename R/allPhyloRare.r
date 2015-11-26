# phylorare: http://davidnipperess.blogspot.co.nz/search/label/Software


if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("tableFile")) stop("table file is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("taxa.group")) taxa.group="assigned"
if(!exists("isPlot")) isPlot = FALSE # by subplot

#source("Modules/init.r")
source("Modules/PhylorareTable.R")
source("Modules/RarefactionPlot.R")

n <- length(matrixNames) 

cat("\nAnalysis: plot the rarefied Phylogenetic Diversity of eDNA community. \n")

env <- getSampleMetaData(isPlot)
env[,"ForestType"] <- gsub(":.*", "", env[,"ForestType"], ignore.case = T)
env[,"ForestType"] <- gsub("x", "unknown", env[,"ForestType"], ignore.case = T)

for (expId in 1:(n-1)) {
  ### eDNA ###
  phylo.rare.df <- getPhyloRareTable(expId, isPlot, rmSingleton, taxa.group)
  
  if (is.null(phylo.rare.df)) {
    cat("\nSkip", taxa.group, "subset from", matrixNames[expId], "because of no tree file.\n")
    next
  }
  
  fname <- paste("phylorare", matrixNames[expId], postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-")
  fileStem <- file.path(workingPath, figDir, fname)
  plotRarefactionTable(phylo.rare.df, env, fileStem)
}
