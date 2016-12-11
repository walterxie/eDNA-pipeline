# Figure 4: Multivariate similarity of non-singleton OTU assemblages

# 16S bacterial (a), 18S protistan (b), 18S fungal (c), 18S animal (d), 26S fungal (e), and COI-300 animal (f)
# by.plot=FALSE to have paired subplots linked
getNMDS <- function(input.names, metric="jaccard", by.plot=FALSE,
                    genes.taxa=list(list("16S","bacteria"),list("18S","protists"),list("18S","fungi"),
                                    list("18S","animals"),list("26S","fungi"),list("ShCO1","animals")) ) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)

  cm.by.subplot.list <- getCommunityList(genes=input.names, genes.taxa=genes.taxa, by.plot=by.plot, 
                              col.ranks=c("superkingdom", "kingdom"), drop.taxa=TRUE )
  env.subplot <- getEnvData(by.plot=by.plot)
  cat("\n")
  
  plot.list <- list()
  for (i in 1:length(cm.by.subplot.list)) {
    cat("NMDS for", names(cm.by.subplot.list)[i], ".\n")
    gg <- ComMA::ggNMDSPlot(t(cm.by.subplot.list[[i]]), env.subplot, colour.id="Elevation", shape.id="Forest.code", link.id="Plot", 
                            distance=metric, shapes=c(15,16,17,0,1,2,5,6,3,4), palette=c("blue", "orange"),
                            shape.levels = c("VS2","VS3","VS5","WF7","WF9","WF11","WF12","WF13","MF20","Unknown"),
                            title = paste0(names(cm.by.subplot.list)[i], ", ", ComMA::simpleCap(metric), "distance"),
                            legend.title.shape="Forest Type", legend.title.colour="Elevation (m)", scale.limits.min=0 )
    plot.list[[names(cm.by.subplot.list)[i]]] <- gg
  }
  
  list( plot.list=plot.list, cm.list=cm.by.subplot.list, metric=metric, by.plot=by.plot  )
}
