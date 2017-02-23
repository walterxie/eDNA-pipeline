# Figure 4: Multivariate similarity of non-singleton OTU assemblages

# 16S bacterial (a), 18S protistan (b), 18S fungal (c), 18S animal (d), 26S fungal (e), and COI-300 animal (f)
# by.plot=FALSE to have paired subplots linked
getNMDS <- function(input.names, metric="jaccard", 
                    genes.taxa=list(list("16S","bacteria"),list("18S","protists"),list("18S","fungi"),
                                    list("18S","animals"),list("26S","fungi"),list("ShCO1","animals")) ) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)

  cm.by.subplot.list <- getCommunityList(genes=input.names, genes.taxa=genes.taxa, by.plot=F )
  cat("\n")
  cm.prep.list <- preprocessCMList(cm.by.subplot.list) 
  cat("\n")
  env.subplot <- getEnvData(by.plot=F)
  cat("\n")
  
  plot.list <- list()
  dist.list <- list()
  require(gg1L)
  for (i in 1:length(cm.prep.list)) {
    cm.name <- names(cm.prep.list)[i]
    cat("NMDS for", cm.name, ".\n")
    cm.prep <- cm.prep.list[[cm.name]]
    
    dsi.dist <- ComMA::getDissimilarity(cm.prep, method=metric)
    dist.list[[cm.name]] <- dsi.dist
    
    gg <- gg1L::ggNMDSPlot(dsi.dist, env.subplot, colour.id="Elevation", link.id="Plot", 
                            shape.id="Forest.code", shapes=c(15,16,17,0,1,2,5,6,3,4), 
                            palette=c("blue", "orange"), text.repel=T, text.size=2.5,
                            shape.levels = c("VS2","VS3","VS5","WF7","WF9","WF11","WF12","WF13","MF20","Unknown"),
                            title = paste0(letters[i], ". ", cm.name, ", ", 
                                           ComMA::simpleCap(metric), " distance"),
                            x.lab="", y.lab="", stress.digits=3, scale.limits.min=0,
                            legend.title.shape="Forest Type", legend.title.colour="Elevation (m)" )
    plot.list[[cm.name]] <- gg
  }
  
  list( plot.list=plot.list, cm.list=cm.prep.list, dist.list=dist.list, metric=metric, by.plot=F  )
}

# besides gt4, give genes.taxa to genes.taxa.list
plotAllNMDS <- function(input.names, genes.taxa.list) {
  require(ggplot2)
  require(gg1L)
  theme_set(theme_bw(base_size=8))
  nmds <- getNMDS(input.names)
  gt4 <- gg1L::grid_arrange_shared_legend(nmds$plot.list, input.list=T, nrow=3, 
                                           legend.position="right", widths=c(0.8, 0.2))
  nmds.list <- list(gt4=gt4)
  
  for (i in 1:length(genes.taxa.list)) {
    gt.name <- names(genes.taxa.list)[i]
    genes.taxa <- genes.taxa.list[[i]]
    nmds <- getNMDS(input.names, genes.taxa=genes.taxa)
    nrow <- round(length(genes.taxa)/2)
    gt <- gg1L::grid_arrange_shared_legend(nmds$plot.list, input.list=T, nrow=nrow, 
                                            legend.position="right", widths=c(0.8, 0.2))
    nmds.list[[gt.name]] <- gt
  }
  return(nmds.list)
}
