


taxa.row.names <- c("ARCHAEA", "BACTERIA", "CHROMISTA", "PROTOZOA", "FUNGI", "PLANTAE", "ANIMALIA", "EUKARYOTA")

getTaxaGroupStatistics <- function(by.plot=TRUE, file.xtable=NULL) {
  source("R/init.R", local=TRUE)
  
  if(!exists("input.names")) stop("input names are missing !")
  output.names <- getOutputNames(input.names)
  
  # data frame for statistics
  tg.otus <- data.frame(stringsAsFactors=FALSE, check.names=FALSE)
  tg.phyla <- data.frame(stringsAsFactors=FALSE, check.names=FALSE)
  
  require(ComMA)
  for (dataId in 1:length(input.names)) {
    # no singletons
    min2=TRUE
    cat("\n", output.names[dataId], "Taxonomic composition,", ifelse(min2, "exclude", "include"), 
        "singletons, samples are based on", ifelse(by.plot, "plot", "subplot"), ".\n") 
    
    cm <- getCommunityMatrix(input.names[dataId], min2=min2, by.plot=by.plot)
    tt <- getTaxaTable(input.names[dataId], taxa.group="assigned")
    
    for (taxaId in 1:length(taxa.row.names)) {
      cat("Taxonomic group", taxa.row.names[taxaId], ".\n") 
      
      tt.sub <- ComMA::subsetTaxaTable(tt, taxa.group=taxa.row.names[taxaId], rank="kingdom")
      
      if (nrow(tt.sub) < 1) {
        tg.otus[taxa.row.names[taxaId],output.names[dataId]] <- 0
        tg.phyla[taxa.row.names[taxaId],output.names[dataId]] <- 0
      } else {
        cm.taxa <- ComMA::mergeCMTaxa(cm, tt.sub, col.ranks = c("kingdom", "phylum"))
        taxa.assign <- ComMA::assignTaxaByRank(cm.taxa)
        
        tg.otus[taxa.row.names[taxaId],output.names[dataId]] <- nrow(cm.taxa)
        tg.phyla[taxa.row.names[taxaId],output.names[dataId]] <- nrow(taxa.assign$phylum)
      }
    }
  }
  
  align.v <- rep("r", ncol(tg.otus) + 1)
  ComMA::printXTable(tg.otus, align = align.v, label = "tab:stats", file = file.xtable,
              caption = paste("") )
  
  align.v <- rep("r", ncol(tg.phyla) + 1)
  ComMA::printXTable(tg.phyla, align = align.v, label = "tab:stats", file = file.xtable,
              caption = paste("") )
  
  list( otus=tg.otus, phyla = tg.phyla )
}
