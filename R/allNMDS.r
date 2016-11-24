# Figure 4: Multivariate similarity of non-singleton OTU assemblages

# 16S bacterial (a), 18S protistan (b), 18S fungal (c), 18S animal (d), 26S fungal (e), and COI-300 animal (f)
# nmds <- getNMDS()
getNMDS <- function(input.names, metric="jaccard",
                    genes.taxa=list(list("16S","bacteria"),list("18S","protists"),list("18S","fungi"),
                                    list("18S","animals"),list("26S","fungi"),list("ShCO1","animals")) ) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)

  cm.list <- getCommunityList(genes=input.names, genes.taxa=genes.taxa, by.plot=F, 
                              col.ranks=c("superkingdom", "kingdom"), drop.taxa=TRUE )
  env.subplot <- getEnvData(by.plot=F)
  cat("\n")
  
  plot.list <- list()
  for (i in 1:length(cm.list)) {
    cat("NMDS for", names(cm.list)[i], ".\n")
    gg <- ComMA::ggNMDSPlot(t(cm.list[[i]]), env.subplot, colour.id="Elevation", shape.id="Forest.code", link.id="Plot", 
                            distance=metric, shapes=c(15,16,17,0,1,2,5,6,3,4), palette=c("blue", "orange"),
                            shape.levels = c("VS2","VS3","VS5","WF7","WF9","WF11","WF12","WF13","MF20","Unknown"),
                            title = paste0(names(cm.list)[i], ", ", ComMA::simpleCap(metric), "distance"),
                            legend.title.shape="Forest Type", legend.title.colour="Elevation (m)", scale.limits.min=0 )
    plot.list[[names(cm.list)[i]]] <- gg
  }
  
  list( plot.list=plot.list, cm.list=cm.list, metric=metric, by.plot=by.plot  )
}
