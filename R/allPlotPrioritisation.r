# Ranking of sample plots by their contributions to the total island biodiversity
# Figure 8: phylogenetic alpha diversity
# Figure S16-S18: Jost biodiversity





prioriPlotByJostDiver <- function(input.names, phylo.tree=NA,
                                  diversities=c("gamma0","gamma1","beta0","beta1","pd.alpha","sp.rich"), 
                                  genes.taxa=list(list("16S","bacteria"),list("18S","animals"),list("18S","fungi"),
                                                  list("18S","protists"), list("26S","animals"),list("26S","fungi"),
                                                  list("26S","protists"),list("ITS","fungi"), list("ShCO1","animals"),
                                                  list("ShCO1","fungi"),list("ShCO1","protists"),
                                                  list("FolCO1","animals"),list("FolCO1","protists"))) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)
  
  cm.list <- getCommunityList(genes=input.names, genes.taxa=genes.taxa, by.plot=T, drop.taxa=TRUE )
  cat("\n")
  
  pd.list <- getPlotPrior(cm.list, input.list=TRUE, is.transposed=FALSE, 
                          phylo.tree=phylo.tree, diversities=diversities)
  return(pd.list)
}








