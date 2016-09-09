# Table 2: Sequence processing and OTU clustering statistics, 
# effective biodiversity and overall taxonomic composition of each amplicon

otus.row.names <- c("Reads", "97\\% OTUs", "Singleton OTUs", "Non-singleton OTUs", "Maximum OTU abundance") 
div.row.names <- c("$^0D_\\alpha$","$^1D_\\alpha$","$^2D_\\alpha$","$^0D_\\beta$","$^1D_\\beta$",
                   "$^2D_\\beta$","$^0D_\\gamma$","$^1D_\\gamma$","$^2D_\\gamma$")

# file.xtable is the file to log xtable
getOTUStatistics <- function(by.plot=TRUE, file.xtable=NULL) {
  source("R/init.R", local=TRUE)
  
  if(!exists("input.names")) stop("input names are missing !")
  output.names <- getOutputNames(input.names)

  # data frame for statistics
  otu.stats <- data.frame(stringsAsFactors=FALSE, check.names=FALSE)
  
  require(ComMA)
  for (dataId in 1:length(input.names)) {
    # inlcude singletons
    min2=FALSE
    cat("\n", output.names[dataId], "OTU clustering statistics,", ifelse(min2, "exclude", "include"), 
        "singletons, samples are based on", ifelse(by.plot, "plot", "subplot"), ".\n") 
    
    cm <- getCommunityMatrix(input.names[dataId], min2=min2, by.plot=by.plot)
    cm.summary <- ComMA::summaryCM.Vector(cm)
    # otus.row.names
    otu.stats[1,output.names[dataId]] <- cm.summary["reads"]
    otu.stats[2,output.names[dataId]] <- cm.summary["OTUs"]
    otu.stats[3,output.names[dataId]] <- cm.summary["singletons"]
    otu.stats[4,output.names[dataId]] <- as.numeric(cm.summary["OTUs"])-as.numeric(cm.summary["singletons"])
    otu.stats[5,output.names[dataId]] <- cm.summary["max.OTU.abundance"]
    row.offset <- 5
    
    # no singletons
    min2=TRUE
    cat("\n", output.names[dataId], "Jost diversities,", ifelse(min2, "exclude", "include"), 
        "singletons, samples are based on", ifelse(by.plot, "plot", "subplot"), ".\n") 
    
    cm <- getCommunityMatrix(input.names[dataId], min2=min2, by.plot=by.plot)
    t.cm <- ComMA::transposeDF(cm)
    # gamma(q=0), alpha(q=0), beta(q=0), gamma(q=1), ..., beta(q=2)
    diversity.v <- ComMA::diversityTable(t.cm, named.vector=T)
    # reorder to match div.row.names
    diversity.v <- diversity.v[c(2,5,8,3,6,9,1,4,7)]
    # div.row.names
    for (i in 1:length(diversity.v)) 
      otu.stats[i+row.offset, output.names[dataId]] <- diversity.v[i]
  }
  rownames(otu.stats) <- c(otus.row.names, div.row.names)
  
  otu.stats <- ComMA::prettyNumbers(otu.stats)
  
  align.v <- rep("r", ncol(otu.stats) + 1)
  ComMA::printXTable(otu.stats, align = align.v, label = "tab:stats", file = file.xtable,
              caption = paste("Sequence processing and OTU clustering statistics,", 
                              "effective biodiversity and overall taxonomic composition of each amplicon") )
  return(otu.stats)
}

