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

pdFilePath <- file.path(workingPath, "data", "pd")
mkdir(pdFilePath) 

cat("\nAnalysis: eDNA community phylogenetic structure. \n")

#env <- getSampleMetaData(isPlot)
#env[,"ForestType"] <- gsub(":.*", "", env[,"ForestType"], ignore.case = T)
#env[,"ForestType"] <- gsub("x", "unknown", env[,"ForestType"], ignore.case = T)

for (expId in 1:(n-1)) {
  ### eDNA ###
  # 16S-assigned-min2, isPlot always FALSE for a valid treeFileStem
  treeFileStem <- tolower(paste(matrixNames[expId], postfix(taxa.group, FALSE, rmSingleton, sep="-"), sep = "-"))
  phyloTree <- getPhyloTree(treeFileStem)
  
  if (is.null(phyloTree)) {
    cat("\nSkip", taxa.group, "subset from", matrixNames[expId], "because of no tree file.\n")
    next
  }
  
  communityMatrix <- getCommunityMatrixT(expId, isPlot, rmSingleton, taxa.group)
  
  pd.alpha <- phylo.alpha(communityMatrix, phyloTree)
  
  pd.phydist <- phylo.phydist(communityMatrix, phyloTree)
  
  pd.mpd <- phylo.mpd(communityMatrix, pd.phydist)
  pd.mntd <- phylo.mntd(communityMatrix, pd.phydist)
  
  pd.beta.dist <- phylo.beta.dist(communityMatrix, pd.phydist)
  
  fileStem <- paste(matrixNames[expId], postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-")
  printResult(pd.alpha, pd.mpd, pd.mntd, pd.beta.dist, pdFilePath, fileStem, tableFile)
    
  fileStem <- paste("pd-beta", matrixNames[expId], postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-")
  plotPDBeta(pd.beta.dist, file.path(workingPath, figDir, fileStem))
}
