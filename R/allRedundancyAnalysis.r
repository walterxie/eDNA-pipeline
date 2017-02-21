# Figure S23 - S27:
# Table S3 - S20:


getEnv <- function(by.plot=FALSE) {
  env <- getEnvData(by.plot=by.plot)
  names(env)[names(env) == 'mean_T_surface'] <- 'Temp.'
  names(env)[names(env) == 'Northness'] <- 'sin.aspect'
  names(env)[names(env) == 'Eastness'] <- 'cos.aspect'
  colnames(env) <- gsub("\\.N", " N", colnames(env))
  colnames(env) <- gsub("\\.P", " P", colnames(env))                                  
  colnames(env) <- gsub("\\.C", " C", colnames(env))                                  
  colnames(env) <- gsub("EC", "E.C.", colnames(env))
  colnames(env) <- gsub("C N.", "C/N ", colnames(env))
  # drop CM30b51 and CM30b58, missing aspect data
  env.prep <- ComMA::preprocessEnv(env, rm.samples=c("CM30b51","CM30b58"), 
                                   sel.env.var=c(4,5,8,9,14:22), log.var=c(5:8,9:11) )
  return(env.prep)
}

# 
getRDA <- function(input.names, by.plot=FALSE, 
                   genes.taxa=list(list("16S","prokaryota"),list("18S","eukaryota"),list("26S","eukaryota"),
                                   list("ITS","eukaryota"),list("ShCO1","eukaryota"),list("FolCO1","eukaryota")) ) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)
  
  cm.list <- getCommunityList(genes=input.names, genes.taxa=genes.taxa, by.plot=by.plot )
  cat("\n")
  # drop CM30b51 and CM30b58, missing aspect data
  cm.prep.list <- preprocessCMList(cm.list, rm.samples=c("CM30b51","CM30b58")) 
  cat("\n")
  env.prep <- getEnvData(by.plot=by.plot)
  cat("\n")
  
  rda.list <- list()
  tcm.list <- list()
  env.list <- list()
  for (i in 1:length(cm.list)) {
    cm.name <- names(cm.list)[i]
    cat("Start RDA analysis for", cm.name, "...\n")
    
    
    
    tcm.env <- ComMA::preprocessRDA(cm.prep.list[[cm.name]], env.prep)
    tcm.list[[cm.name]] <- tcm.env$tcm
    env.list[[cm.name]] <- tcm.env$env
    
    rda <- ComMA::proceedRDA(tcm.env$tcm, tcm.env$env)
    rda.list[[cm.name]] <- rda
  }
  
  list(rda.list=rda.list, tcm.list=tcm.list, env.list=env.list)
}



