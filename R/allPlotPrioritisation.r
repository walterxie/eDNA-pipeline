# Ranking of sample plots by their contributions to the total island biodiversity
# Figure 8: phylogenetic alpha diversity
# Figure S16-S18: Jost biodiversity

prioriPlotByDiversities <- function(input.names, 
                                    diversities=c("gamma0","gamma1","beta0","beta1","pd.alpha","sp.rich"), 
                                    genes.taxa=list(list("16S","prokaryota"),list("18S","eukaryota"),list("26S","eukaryota"),
                                                    list("ITS","eukaryota"),list("ShCO1","eukaryota"),list("FolCO1","eukaryota")) ) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)
  
  cm.by.plot.list <- getCommunityList(genes=input.names, genes.taxa=genes.taxa, by.plot=T, drop.taxa=TRUE )
  cat("\n")
  # need for "pd.alpha","sp.rich"
  tre.list <- getTreeList(genes=input.names, genes.taxa=genes.taxa)
  cat("\n")
  
  cm.list <- list()
  for (i in 1:length(cm.by.plot.list)) {
    cm.name <- names(cm.by.plot.list)[i]
    cat("Preprocess ", cm.name, ".\n")
    cm.prep <- ComMA::preprocessCM(cm.by.plot.list[[cm.name]], min.abund=5, mean.abund.thr=0.025 )
    cm.list[[cm.name]] <- cm.prep
  }
  cat("\n")
  
  plot.prior.list <- ComMA::getPlotPrior(cm.list, is.transposed=FALSE, tre.list=tre.list, taxa.match=FALSE,
                                         diversities=diversities)
  #plot.prior.list <- ComMA::getPlotPrior(cm.list, is.transposed=FALSE, diversities=diversities)
  cat("\n")
  pp.df.list <- ComMA::mergePlotPriorListOfDF(plot.prior.list)
  return(pp.df.list)
}

# give ranks to pp.df.list, such as pp.df.list[["rank"]][["pd.alpha"]]
plotAllHeatmaps <- function(ranks.list, env.plot, pattern="\\..*", replacement="", y.lab="Sample plot", 
                            x.lab="Amplicon dataset", grid.widths = c(8, 2)) {
  theme_set(theme_bw(base_size=8))
  hm.list <- list()
  for (i in 1:length(ranks.list)) {
    gt.name <- names(ranks.list)[i]
    ranks <- replaceColNames(ranks.list[[i]], pattern=pattern, replacement=replacement)
    gt <- ComMA::plotPrioritisation.Phyloseq(ranks, env.plot, x.lab=x.lab, y.lab=y.lab, grid.widths=grid.widths)
    hm.list[[gt.name]] <- gt$heatmap
  }
  return(hm.list)
}





