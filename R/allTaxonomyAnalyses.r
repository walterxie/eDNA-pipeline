# Figure 2 The number of sequences and OTUs detected for each phylum and amplicon. 
# Diamonds represent the number of sequences, 
# open circles the number of OTUs including singleton OTUs, 
# and filled circles the number of OTUs excluding singleton OTUs.

# all.counts.sums <- getAllCountsSums()
getAllCountsSums <- function(by.plot=TRUE, file.xtable=NULL, invalid.char=FALSE, init=TRUE) {
  if (init) source("R/init.R", local=TRUE)
  
  if(!exists("input.names")) stop("input names are missing !")
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
    tt$kingdom <- gsub("CHROMISTA|PROTOZOA", "PROTISTS", tt$kingdom)
    
    # preprocess taxonomy
    cm.taxa <- ComMA::mergeCMTaxa(cm, tt, col.ranks = c("superkingdom", "kingdom", "phylum"))
    cm.taxa.list[[output.names[data.id]]] <- cm.taxa
  }
  cat("\n")
  
  all.counts.sums <- ComMA::sumReadsOTUs(cm.taxa.list, taxa.rank="phylum", group.rank="kingdom", 
                                         col.ranks=c("superkingdom", "kingdom", "phylum"), 
                                         pdf.width = 340, pdf.height = 260)

  return(all.counts.sums)
}

