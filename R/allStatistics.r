# Table 2: Sequence processing and OTU clustering statistics, 
# effective biodiversity and overall taxonomic composition of each amplicon

# file.xtable is the file to log xtable
getOTUStatistics <- function(by.plot=TRUE, file.xtable=NULL, invalid.char=FALSE) {
  source("R/init.R", local=TRUE)
  
  if(!exists("input.names")) stop("input names are missing !")
  output.names <- getOutputNames(input.names)

  cm.list <- list()
  cm.min2.list <- list()
  require(ComMA)
  for (data.id in 1:length(input.names)) {
    # inlcude singletons
    min2=FALSE
    cat("\n", output.names[data.id], "OTU clustering statistics,", ifelse(min2, "exclude", "include"), 
        "singletons, samples are based on", ifelse(by.plot, "plot", "subplot"), ".\n") 
    cm <- getCommunityMatrix(input.names[data.id], min2=min2, by.plot=by.plot)
    cm.list[[output.names[data.id]]] <- cm
    
    # no singletons
    min2=TRUE
    cat("\n", output.names[data.id], "Jost diversities,", ifelse(min2, "exclude", "include"), 
        "singletons, samples are based on", ifelse(by.plot, "plot", "subplot"), ".\n") 
    cm <- getCommunityMatrix(input.names[data.id], min2=min2, by.plot=by.plot)
    cm.min2.list[[output.names[data.id]]] <- cm
  }
  cat("\n")
  
  # 1. OTU summary
  otu.stats <- ComMA::summaryOTUs(cm.list, input.list=T)
  # hard code to keep rows we need
  otu.stats <- otu.stats[c(1,2,4,10,6),]
  rownames(otu.stats) <- c("Reads", "97\\% OTUs", "Singleton OTUs", 
                           "Non-singleton OTUs", "Maximum OTU abundance")
  # 2. Jost diversity summary
  div.stats <- ComMA::summaryDiversity(cm.min2.list, input.list=T, row.order=c(2,5,8,3,6,9,1,4,7))
  otu.stats <- rbind(otu.stats, div.stats)
  
  # 3. xtable
  align.v <- rep("r", ncol(otu.stats) + 1)
  ComMA::printXTable(otu.stats, align = align.v, label = "tab:stats", file = file.xtable, invalid.char=invalid.char,
              caption = paste("Sequence processing and OTU clustering statistics,", 
                              "effective biodiversity and overall taxonomic composition of each amplicon") )
  return(otu.stats)
}

