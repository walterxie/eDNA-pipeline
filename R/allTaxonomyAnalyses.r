# Figure 2 The number of sequences and OTUs detected for each phylum and amplicon. 
# Diamonds represent the number of sequences, 
# open circles the number of OTUs including singleton OTUs, 
# and filled circles the number of OTUs excluding singleton OTUs.

# all.counts.sums <- getAllCountsSums()
getAllCountsSums <- function(input.names, by.plot=TRUE) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)
  output.names <- getOutputNames(input.names)
  
  cm.taxa.list <- list()
  require(ComMA)
  for (data.id in 1:length(input.names)) {
    # inlcude singletons
    min2=FALSE
    cat("\n", output.names[data.id], "OTU clustering statistics,", ifelse(min2, "exclude", "include"), 
        "singletons, samples are based on", ifelse(by.plot, "plot", "subplot"), ".\n") 
    cm <- getCommunityMatrix(input.names[data.id], min2=min2, by.plot=by.plot)
    tt <- getTaxaTable(input.names[data.id], taxa.group="all")
    # adjust taxonomy?
    #tt$kingdom <- gsub("CHROMISTA|PROTOZOA", "PROTISTS", tt$kingdom)
    
    # no preprocess to keep original names
    cm.taxa <- ComMA::mergeCMTaxa(cm, tt, preprocess=F, col.ranks = c("superkingdom", "kingdom", "phylum"))
    cm.taxa.list[[output.names[data.id]]] <- cm.taxa
  }
  cat("\n")
  
  # include figure and table
  all.counts.sums <- ComMA::summReadsOTUsPipeline(cm.taxa.list, taxa.rank="phylum", group.rank="kingdom", 
                                           col.ranks=c("superkingdom", "kingdom", "phylum"), 
                                           gene.levels=c("16S", "18S", "26S", "ITS", "COI-300", "COI-650"),
                                           group.levels=c("ARCHAEA","BACTERIA","EUKARYOTA","PROTOZOA","CHROMISTA",
                                                          "FUNGI","PLANTAE","ANIMALIA","Unknown"),
                                           palette=c("orange","red","blue","steelblue","skyblue",
                                                     "purple","green2","green4","grey"),
                                           legend.title = "(Super)Kingdom")
  
  return(all.counts.sums)
}

printAllCountsSums <- function(all.counts.sums, taxa.rank="phylum", file.xtable=NULL, invalid.char=FALSE) {
  align.v <- rep("r", ncol(all.counts.sums$all.counts.sums) + 1)
  ComMA::printXTable(all.counts.sums$all.counts.sums, align = align.v, label = "tab:counts:sums", 
                     file = file.xtable, invalid.char=invalid.char,
                     caption = paste("The summary of number of reads and OTUs assigned to", taxa.rank) )
}
