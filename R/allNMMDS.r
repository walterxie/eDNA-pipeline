# Figure 4: Multivariate similarity of non-singleton OTU assemblages

# 16S bacterial (a), 18S protistan (b), 18S fungal (c), 18S animal (d), 26S fungal (e), and COI-300 animal (f)
# nmds <- getNMDS()




# nmds <- getNMDS()
getNMDS <- function(input.names, by.plot=TRUE, file.xtable=NULL, invalid.char=FALSE) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)
  output.names <- getOutputNames(input.names)
  
  taxa.group <- c("ARCHAEA", "BACTERIA", "CHROMISTA", "PROTOZOA", "FUNGI", "PLANTAE", "ANIMALIA", "EUKARYOTA")
  
  # data frame for statistics
  cm.taxa.list <- list()
  require(ComMA)
  for (data.id in 1:length(input.names)) {
    # no singletons
    min2=TRUE
    cat("\n", output.names[data.id], "Taxonomic composition,", ifelse(min2, "exclude", "include"), 
        "singletons, samples are based on", ifelse(by.plot, "plot", "subplot"), ".\n") 
    
    cm <- getCommunityMatrix(input.names[data.id], min2=min2, by.plot=by.plot)
    tt <- getTaxaTable(input.names[data.id], taxa.group="assigned")
    
    cm.taxa <- ComMA::mergeCMTaxa(cm, tt, col.ranks = c("superkingdom", "kingdom", "phylum"))
    cm.taxa.list[[output.names[data.id]]] <- cm.taxa
  }
  cat("\n")
  
 
  
  list( otus=ta.gr.stats$otus, merged=merged.stats, phyla=ta.gr.stats$rank.count, 
        phyla.list=ta.gr.stats$count.rank.df.list  )
}
