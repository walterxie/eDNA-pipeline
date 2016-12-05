# Figure S23 - S27:
# Table S3 - S20:

# env.subplot <- getEnvData(by.plot=F)
getRDA <- function(input.names, env.subplot, metric="jaccard", 
                   genes.taxa=list(list("16S","prokaryotes"),list("18S","eukaryotes"),list("26S","eukaryotes"),
                                   list("ITS","eukaryotes"),list("ShCO1","eukaryotes"),list("FolCO1","eukaryotes"))) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)
  
  cm.list <- getCommunityList(genes=input.names, genes.taxa=genes.taxa, by.plot=F, drop.taxa=TRUE )
  cat("\n")
  
  
  
  
}



