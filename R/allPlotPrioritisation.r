# Ranking of sample plots by their contributions to the total island biodiversity
# Figure 8: phylogenetic alpha diversity
# Figure S16-S18: Jost biodiversity

prioriPlotByDiversities <- function(input.names, 
                                    diversities=c("gamma0","gamma1","beta0","beta1","pd.alpha","sp.rich"), 
                                    genes.taxa=list(list("16S","prokaryota"),list("18S","eukaryota"),list("26S","eukaryota"),
                                                    list("ITS","eukaryota"),list("ShCO1","eukaryota"),list("FolCO1","eukaryota")) ) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)
  
  cm.list <- getCommunityList(genes=input.names, genes.taxa=genes.taxa, by.plot=T, drop.taxa=TRUE )
  cat("\n")
  # need for "pd.alpha","sp.rich"
  tre.list <- getTreeList(genes=input.names, genes.taxa=genes.taxa)
  cat("\n")
  
  plot.prior.list <- ComMA::getPlotPrior(cm.list, is.transposed=FALSE, tre.list=tre.list, taxa.match=FALSE,
                                         diversities=diversities)
  #plot.prior.list <- ComMA::getPlotPrior(cm.list, is.transposed=FALSE, diversities=diversities)
  cat("\n")
  pp.df.list <- ComMA::mergePlotPriorListOfDF(plot.prior.list)
  return(pp.df.list)
}








