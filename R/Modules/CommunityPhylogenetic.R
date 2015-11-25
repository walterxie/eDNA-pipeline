# Picante: R tools for integrating phylogenies and ecology

library(picante)
library(cluster)

source("Modules/init.R", local=TRUE)

# Community phylogenetic structure
# communityMatrix: rows are simples
# phyloTree: rooted tree of phylo object
# cm.env: rows are simples, and must be same as rownames(communityMatrix) inlcuding order
# treepathFileStem: tree file name but not include .tre
# tableFile: latex file
# verbose default TRUE
phylo.alpha <- function(communityMatrix, phyloTree, ORD.RES=function(res) {res[order(rownames(res)),]}, verbose=TRUE) {
  if ( ! all( sort(colnames(communityMatrix)) == sort(phyloTree$tip.label) ) ) 
    stop( paste("Community OTU names do not match tree tip labels") )
  
  if(verbose) {
    cat("Input community", nrow(communityMatrix), "samples", ncol(communityMatrix), "OTUs", 
        ", phylogenetic tree with", length(phyloTree$tip.label), "tips and", phyloTree$Nnode, "internal nodes.\n") 
    cat("Analysis: Faith's phylogenetic alpha diversity.\n") 
  }
  
  # phylogenetic alpha diversity (PD) index proposed by Faith (1992)
  pd.result <- pd(communityMatrix, phyloTree, include.root = TRUE)
  
  pd.result <- ORD.RES(pd.result)
  pd.result
}
  
phylo.phydist <- function(communityMatrix, phyloTree) {
  # cophenetic distances for a hierarchical clustering
  phydist <- cophenetic(phyloTree)
  phydist
}

phylo.mpd <- function(communityMatrix, phydist, ORD.RES=function(res) {res[order(rownames(res)),]}, verbose=TRUE) {
  if(verbose) 
    cat("Analysis: mean pairwise distance (MPD).\n")
  
  # When used with a phylogenetic distance matrix, equivalent to -1 times the Nearest Taxon Index (NTI).
  # MPD: standardized effect size of mean pairwise distances in communities
  ses.mpd.result <- ses.mpd(communityMatrix, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 99)
  ses.mpd.result <- ORD.RES(ses.mpd.result)
  ses.mpd.result
}

phylo.mntd <- function(communityMatrix, phydist, ORD.RES=function(res) {res[order(rownames(res)),]}, verbose=TRUE) {
  if(verbose) 
    cat("Analysis: mean nearest taxon distance (MNTD).\n")

  # MNTD: standardized effect size of mean nearest taxon distances in communities
  ses.mntd.result <- ses.mntd(communityMatrix, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 99)
  ses.mntd.result <- ORD.RES(ses.mntd.result)
  ses.mntd.result
}

phylo.beta.dist <- function(communityMatrix, phydist) {
  if(verbose) 
    cat("Analysis: MPD phylogenetic beta diversity.\n")
  
  # MPD (mean pairwise distance) separating taxa in two communities, phylogenetic beta diversity (Steven Kembel)
  comdist.result <- comdist(communityMatrix, phydist)
  comdist.result
}

# create file for intermediate data
# pathFileStem = filePath + fileStem
printResult <- function(pd.alpha, pd.mpd, pd.mntd, pd.beta.dist, filePath, fileStem, tableFile=NULL) {
  fname <- file.path(filePath, paste(fileStem, "pd-alpha", sep = "-"))
  write.cm(pd.alpha, paste(fname, "csv", sep = "."))
  printXTable(pd.alpha, 
              caption = paste("Faith's phylogenetic alpha diversity of", fileStem), 
              label = paste("tab:pd:alpha", fileStem, sep = ":"), file=tableFile)    
  
  
  fname <- file.path(filePath, paste(fileStem, "mpd", sep = "-"))
  write.cm(pd.mpd, paste(fname, "csv", sep = "."))
  printXTable(pd.mpd, 
              caption = paste("The mean pairwise distance (MPD) between all species in each community of", fileStem), 
              label = paste("tab:mpd:", fileStem, sep = ":"), file=tableFile)    
  
  fname <- file.path(filePath, paste(fileStem, "mntd", sep = "-"))
  write.cm(pd.mntd, paste(fname, "csv", sep = "."))
  printXTable(pd.mntd, 
              caption = paste("The mean nearest taxon distance (MNTD)", 
                              "separating each species in the community from its closest relative of", fileStem), 
              label = paste("tab:mntd:", fileStem, sep = ":"), file=tableFile) 
  
  comdist.m <- as.matrix(pd.beta.dist)
  
  fname <- file.path(filePath, paste(fileStem, "pd-beta", sep = "-"))
  write.cm(comdist.m, paste(fname, "csv", sep = "."))
  
  comdist.m <- comdist.m[order(rownames(comdist.m)),]
  comdist.m <- comdist.m[,order(colnames(comdist.m))]
  comdist.m[upper.tri(comdist.m)] <- NA
  comdist.m[comdist.m=="0"] <- NA
  comdist.m <- comdist.m[-1,-ncol(comdist.m)]
  
  printXTable(comdist.m, 
              caption = paste("Phylogenetic beta diversity of", fileStem), 
              label = paste("tab:pd:beta", fileStem, sep = ":"), file=tableFile)  
}
  
plotPDBeta <- function(pd.beta.dist, pathFileStem) {
  comdist.clusters <- hclust(pd.beta.dist)
  
  pdf(paste(pathFileStem, "pdf", sep = "."), width=10, height=5)
  plot(comdist.clusters, xlab="", sub ="", main=paste("Phylogenetic beta diversity"))
  invisible(dev.off())
}


