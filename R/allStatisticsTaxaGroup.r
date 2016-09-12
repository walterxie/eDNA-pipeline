# Table 2: Sequence processing and OTU clustering statistics, 
# effective biodiversity and overall taxonomic composition of each amplicon

getTaxaGroupStatistics <- function(by.plot=TRUE, file.xtable=NULL) {
  source("R/init.R", local=TRUE)
  
  if(!exists("input.names")) stop("input names are missing !")
  output.names <- getOutputNames(input.names)
  
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
    
    cm.taxa <- ComMA::mergeCMTaxa(cm, tt, col.ranks = c("kingdom", "phylum"))
    cm.taxa.list[[output.names[data.id]]] <- cm.taxa
  }
  
  cat("\n")
  taxa.group <- c("ARCHAEA", "BACTERIA", "CHROMISTA", "PROTOZOA", "FUNGI", "PLANTAE", "ANIMALIA", "EUKARYOTA")
  # remove every rows containing "unclassified".
  ta.gr.stats <- ComMA::summaryTaxaGroup(cm.taxa.list, input.list=T, unclassified=3, taxa.group=taxa.group,
                                         group.rank="kingdom", count.rank="phylum")
  
  align.v <- rep("r", ncol(ta.gr.stats$otus) + 1)
  ComMA::printXTable(ta.gr.stats$otus, align = align.v, label = "tab:stats", file = file.xtable,
              caption = paste("") )
  
  align.v <- rep("r", ncol(ta.gr.stats$rank.count) + 1)
  ComMA::printXTable(ta.gr.stats$rank.count, align = align.v, label = "tab:stats", file = file.xtable,
              caption = paste("") )
  
  list( otus=ta.gr.stats$otus, phyla = ta.gr.stats$rank.count )
}
