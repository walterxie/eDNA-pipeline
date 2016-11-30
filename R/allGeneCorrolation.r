


# mantel <- getMantel()
getMantelAndProcrustes <- function(input.names, metric="jaccard",
                    genes.taxa=list(list("16S","all"),list("18S","all"),list("26S","all"),
                                    list("ITS","all"),list("ShCO1","all"),list("FolCO1","all")) ) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)
  
  cm.list <- getCommunityList(genes=input.names, genes.taxa=genes.taxa, by.plot=F, 
                              col.ranks=c("superkingdom", "kingdom"), drop.taxa=TRUE )
  cat("\n")
  
  dissim <- ComMA::getDissimilarityList(cm.list, metric=metric)
  mantel <- ComMA::mantelComparison(dissim$dist.list)
  procrustes <- ComMA::procrustesComparison(dissim$dist.list)
  
  mantel.tri <- ComMA::getTriMatrix(mantel$m.df) # Mantel stats
  prot.tri <- ComMA::getTriMatrix(procrustes$prot.df) # Procrustes stats
  corrs <- ComMA::combineTriMatrix(mantel.tri, prot.tri)
  
  list( corrs=corrs, mantel.tri=mantel.tri, prot.tri=prot.tri,
        mantel=mantel, procrustes=procrustes, dist.list=dissim$dist.list, 
        metric=metric, msg="Mantel in upper triangle, Procrustes in lower"  )
}
#ggNMDSPlot(dist(mantel.tri), text.or.point=1, text.size=5, title="") + expand_limits(x = c(-1, 0.8), y=c(-0.8, 0.8))
#ggNMDSPlot(dist(prot.tri), text.or.point=1, text.size=5, title="") + expand_limits(x = c(-0.6, 1), y=c(-0.8, 0.8))

