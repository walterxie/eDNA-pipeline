# Picante: R tools for integrating phylogenies and ecology

library(picante)
library(cluster)

source("Modules/init.R", local=TRUE)

# Community phylogenetic structure
# communityMatrix: rows are simples
# phyloTree: rooted tree of phylo object
# cm.env: rows are simples, and must be same as rownames(communityMatrix) inlcuding order
# treeFileStem: tree file name but not include .tre
# tableFile: latex file
# verbose default TRUE
comm.phylo.struc <- function(communityMatrix, phyloTree, treeFileStem, tableFile=NULL, verbose=TRUE, 
                             ORD.RES=function(res) {res[order(rownames(res)),]}) {
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
  printXTable(pd.result, 
              caption = paste("Faith's phylogenetic alpha diversity of", treeFileStem), 
              label = paste("tab:pd:alpha", treeFileStem, sep = ":"), file=tableFile)    
  
  if(verbose) 
    cat("Analysis: MPD and MNTD.\n")
  
  # cophenetic distances for a hierarchical clustering
  phydist <- cophenetic(phyloTree)
  
  # When used with a phylogenetic distance matrix, equivalent to -1 times the Nearest Taxon Index (NTI).
  # MPD: standardized effect size of mean pairwise distances in communities
  ses.mpd.result <- ses.mpd(communityMatrix, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 99)

  ses.mpd.result <- ORD.RES(ses.mpd.result)
  printXTable(ses.mpd.result, 
              caption = paste("The mean pairwise distance (MPD) between all species in each community of", treeFileStem), 
              label = paste("tab:mpd:", treeFileStem, sep = ":"), file=tableFile)    
  
  # MNTD: standardized effect size of mean nearest taxon distances in communities
  ses.mntd.result <- ses.mntd(communityMatrix, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 99)
  
  ses.mntd.result <- ORD.RES(ses.mntd.result)
  printXTable(ses.mntd.result, 
              caption = paste("The mean nearest taxon distance (MNTD)", 
                              "separating each species in the community from its closest relative of", treeFileStem), 
              label = paste("tab:mntd:", treeFileStem, sep = ":"), file=tableFile)    
  
  if(verbose) 
    cat("Analysis: MPD phylogenetic beta diversity.\n")
  
  # MPD (mean pairwise distance) separating taxa in two communities, phylogenetic beta diversity (Steven Kembel)
  comdist.result <- comdist(communityMatrix, phydist)
  
  comdist.m <- as.matrix(comdist.result)
  comdist.m <- comdist.m[order(rownames(comdist.m)),]
  comdist.m <- comdist.m[,order(colnames(comdist.m))]
  comdist.m[upper.tri(comdist.m)] <- NA
  comdist.m[comdist.m=="0"] <- NA
  comdist.m <- comdist.m[-1,-ncol(comdist.m)]
  
  printXTable(comdist.m, 
              caption = paste("Phylogenetic beta diversity of", treeFileStem), 
              label = paste("tab:pd:beta", treeFileStem, sep = ":"), file=tableFile)    
  
  comdist.clusters <- hclust(comdist.result)
  
  fname <- paste("pd-beta-", treeFileStem, ".pdf", sep = "")
  pdf(file.path(workingPath, figDir, fname), width=10, height=5)
  plot(comdist.clusters, xlab="", sub ="", main=paste("Phylogenetic beta diversity"))
  invisible(dev.off())
}


