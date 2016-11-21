# Table 2: Sequence processing and OTU clustering statistics, 
# effective biodiversity and overall taxonomic composition of each amplicon

# The (EUKARYOTA) category is the number of OTUs that are identified as eukaryotes but aren't identified to kingdom level.
# use phyla.list to check if the assigned phyla in each group are correct
# tg.stats <- getTaxaGroupStatistics()
getTaxaGroupStatistics <- function(input.names, by.plot=TRUE, file.xtable=NULL, invalid.char=FALSE) {
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
  
  # No EUKARYOTA, remove every rows containing "unclassified".
  ta.gr.stats <- ComMA::summaryTaxaGroup(cm.taxa.list, input.list=T, unclassified=3, 
                                         taxa.group=taxa.group[-length(taxa.group)],
                                         group.rank="kingdom", count.rank="phylum")
  merged.stats <- ComMA::combineTwoDF(ta.gr.stats$otus, ta.gr.stats$rank.count)
  
  # OTUs only identified to EUKARYOTA
  ta.gr.stats2 <- ComMA::summaryTaxaGroup(cm.taxa.list, input.list=T, unclassified=3, pretty.numbers=F,
                                         taxa.group=taxa.group[length(taxa.group)],
                                         group.rank="superkingdom", count.rank="kingdom")
  only.euk <- ta.gr.stats2$otus - sapply(ta.gr.stats2$count.rank.df.list$OTUs, sum)
  merged.stats <- rbind(merged.stats, ComMA::prettyNumbers(only.euk, digits = 0))
  
  cat("\n")
  align.v <- rep("r", ncol(merged.stats) + 1)
  ComMA::printXTable(merged.stats, align = align.v, label = "tab:tgroup:stats", 
                     file = file.xtable, invalid.char=invalid.char,
              caption = paste("Overall taxonomic composition of each amplicon") )
  
  list( otus=ta.gr.stats$otus, merged=merged.stats, phyla=ta.gr.stats$rank.count, 
        phyla.list=ta.gr.stats$count.rank.df.list  )
}
