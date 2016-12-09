# Figure S23 - S27:
# Table S3 - S20:

# 
getRDA <- function(input.names, by.plot=FALSE, metric="jaccard", 
                   genes.taxa=list(list("16S","prokaryotes"),list("18S","eukaryotes"),list("26S","eukaryotes"),
                                   list("ITS","eukaryotes"),list("ShCO1","eukaryotes"),list("FolCO1","eukaryotes"))) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)
  
  cm.list <- getCommunityList(genes=input.names, genes.taxa=genes.taxa, by.plot=by.plot, drop.taxa=TRUE )
  cat("\n")
  env <- getEnvData(by.plot=by.plot)
  names(env)[names(env) == 'mean_T_surface'] <- 'Temp.'
  names(env)[names(env) == 'Northness'] <- 'sin.aspect'
  names(env)[names(env) == 'Eastness'] <- 'cos.aspect'
  colnames(env) <- gsub("\\.N", " N", colnames(env))
  colnames(env) <- gsub("\\.P", " P", colnames(env))                                  
  colnames(env) <- gsub("\\.C", " C", colnames(env))                                  
  colnames(env) <- gsub("EC", "E.C.", colnames(env))
  colnames(env) <- gsub("C N.", "C/N ", colnames(env))
  
  rda.list <- list()
  for (i in 1:length(cm.list)) {
    cm <- cm.list[[i]]
    cat("Start RDA analysis for", names(cm.list)[i], "with", ncol(cm), "samples initially.\n")
    
    # drop CM30b51 and CM30b58, missing aspect data
    tcm.env <- ComMA::preprocessRDA(cm, env, is.transposed=FALSE, env.var=c(4,5,8,9,14:22), rm="CM30b51|CM30b58")
    rda <- ComMA::proceedRDA(tcm.env$tcm, tcm.env$env)
    
    
  }
  
  list(rda.list=rda.list, cm.list=cm.list, env=env)
}



