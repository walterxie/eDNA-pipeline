# Figure S23 - S27:
# Table S3 - S20:

# 
getRDA <- function(input.names, by.plot=FALSE, 
                   genes.taxa=list(list("16S","prokaryota"),list("18S","eukaryota"),list("26S","eukaryota"),
                                   list("ITS","eukaryota"),list("ShCO1","eukaryota"),list("FolCO1","eukaryota")) ) {
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
  tcm.list <- list()
  env.list <- list()
  for (i in 1:length(cm.list)) {
    cm <- cm.list[[i]]
    cat("Start RDA analysis for", names(cm.list)[i], "with", ncol(cm), "samples initially.\n")
    
    # drop CM30b51 and CM30b58, missing aspect data
    tcm.env <- ComMA::preprocessRDA(cm, env, is.transposed=FALSE, rm="CM30b51|CM30b58", 
                                    env.var=c(4,5,8,9,14:22), log.var=c(5:8,9:11))
    tcm.list[[names(cm.list)[i]]] <- tcm.env$tcm
    env.list[[names(cm.list)[i]]] <- tcm.env$env
    
    rda <- ComMA::proceedRDA(tcm.env$tcm, tcm.env$env)
    rda.list[[names(cm.list)[i]]] <- rda
  }
  
  list(rda.list=rda.list, tcm.list=tcm.list, env.list=env.list)
}



