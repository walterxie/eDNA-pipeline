# Table 2: Sequence processing and OTU clustering statistics, 
# effective biodiversity and overall taxonomic composition of each amplicon

# use phyla.list to check if the assigned phyla in each group are correct
getTaxaGroupStatistics <- function(by.plot=TRUE, file.xtable=NULL, invalid.char=FALSE) {
  source("R/init.R", local=TRUE)
  taxa.group <- c("ARCHAEA", "BACTERIA", "CHROMISTA", "PROTOZOA", "FUNGI", "PLANTAE", "ANIMALIA", "EUKARYOTA")
  
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
    
    cm.taxa <- ComMA::mergeCMTaxa(cm, tt, col.ranks = c("superkingdom", "kingdom", "phylum"))
    cm.taxa.list[[output.names[data.id]]] <- cm.taxa
  }
  
  cat("\n")
  # remove every rows containing "unclassified".
  ta.gr.stats <- ComMA::summaryTaxaGroup(cm.taxa.list, input.list=T, unclassified=3, 
                                         taxa.group=taxa.group[1:length(taxa.group)-1],
                                         group.rank="kingdom", count.rank="phylum")
  # add superkingdom
  ta.gr.stats2 <- ComMA::summaryTaxaGroup(cm.taxa.list, input.list=T, unclassified=3, 
                                         taxa.group=taxa.group[length(taxa.group)],
                                         group.rank="superkingdom", count.rank="phylum")
  otus <- rbind(ta.gr.stats$otus, ta.gr.stats2$otus)
  phyla <-rbind(ta.gr.stats$rank.count, ta.gr.stats2$rank.count)
  merged.stats <- ComMA::merge2DF(otus, phyla)
  
  # append phyla.list
  phyla.list <- append(ta.gr.stats$count.rank.df.list, ta.gr.stats2$count.rank.df.list)
  
  align.v <- rep("r", ncol(merged.stats) + 1)
  ComMA::printXTable(merged.stats, align = align.v, label = "tab:tgroup:stats", 
                     file = file.xtable, invalid.char=invalid.char,
              caption = paste("Overall taxonomic composition of each amplicon") )
  
  list( otus=otus, merged=merged.stats, phyla=phyla,phyla.list=phyla.list  )
}
