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
                                   log.var=c(14:20), sel.env.var=c(4,5,8,9,14:22) )
  return(env.prep)
}

# rda.list
getRDAList <- function(input.names, by.plot=FALSE, 
                   genes.taxa=list(list("16S","prokaryota"),list("18S","eukaryota"),list("26S","eukaryota"),
                                   list("ITS","eukaryota"),list("ShCO1","eukaryota"),list("FolCO1","eukaryota")) ) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)
  
  cm.list <- getCommunityList(genes=input.names, genes.taxa=genes.taxa, by.plot=by.plot )
  cat("\n")
  # drop CM30b51 and CM30b58, missing aspect data
  cm.prep.list <- preprocessCMList(cm.list, rm.samples=c("CM30b51","CM30b58"), min.abund=5, mean.abund.thr=0.025) 
  cat("\n")
  env.prep <- getEnv(by.plot=by.plot)
  sapply(env.prep, class)
  env.prep <- ComMA::convertType(env.prep)
  cat("\n")
  
  rda.list <- list()
  tcm.list <- list()
  env.list <- list()
  for (i in 1:length(cm.prep.list)) {
    cm.name <- names(cm.prep.list)[i]
    cat("\nStart RDA analysis for", cm.name, "...\n")
    
    tcm.env <- ComMA::matchCMEnv(cm.prep.list[[cm.name]], env.prep)
    tcm.list[[cm.name]] <- tcm.env$tcm
    env.list[[cm.name]] <- tcm.env$env
    
    rda <- ComMA::doRDA(tcm.env$tcm, tcm.env$env)
    rda.list[[cm.name]] <- rda
  }
  
  list(rda.list=rda.list, tcm.list=tcm.list, env.list=env.list)
}

# "pick" to choose the model from "reduced", "forward", "backward"
plotAllRDA <- function(rda.list, env.list, pick="backward") {
  require(ggplot2)
  require(gg1L)
  theme_set(theme_bw(base_size=8))
  
  rda.pl.list <- list()
  for (i in 1:length(rda.list)) {
    cm.name <- names(rda.list)[i]
    rda.pl <- plotRDA(rda.list[[i]][[pick]], env.list[[i]], scale.limits.min=0,
                      title=paste0(letters[i], ". ", cm.name, ", ",  "Jaccard distance"))
    rda.pl.list[[cm.name]] <- rda.pl$plot
  }
  nrow <- round(length(rda.list)/2)
  gt <- gg1L::grid_arrange_shared_legend(rda.pl.list, input.list=T, nrow=nrow, 
                                          legend.position="right", widths=c(0.8, 0.2))
  list(gt=gt, rda.pl.list=rda.pl.list, pick=pick)
}

printAllRDA <- function(rda.list, file.xtable=NULL, invalid.char=FALSE) {
  for (i in 1:length(rda.list)) {
    cm.name <- names(rda.list)[i]
    gene.taxa <- unlist(strsplit(cm.name, split=" "))
    ComMA::printXTable.RDA(rda.list[[i]], matrix.name=gene.taxa[1], taxa.group=gene.taxa[2], 
                           table.file=file.xtable, invalid.char=invalid.char)
  }
}



