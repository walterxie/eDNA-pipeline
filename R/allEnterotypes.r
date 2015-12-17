library(gplots)
library(ggplot2)
library(grid)
library(RColorBrewer)


if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("isPlot")) isPlot = FALSE # by subplot
if(!exists("taxa.group")) taxa.group="all"
if(!exists("rankLevel")) rankLevel="phylum"
if(!exists("groupLevel")) groupLevel="kingdom"

n <- length(matrixNames) 

source("Modules/init.r")
source("Modules/Enterotypes.R")

######## eDNA by plot/subplot #########
for (expId in 1:(n-1)) {
  taxaAssgReads <- getTaxaAssgReads(expId, isPlot, rmSingleton, rankLevel, groupLevel, taxa.group) 
  
  if (nrow(taxaAssgReads) < 1) {
    cat("\nSkip", taxa.group, "subset taxonomic summary from", matrixNames[expId], ".\n") 
    next
  }
  
  taxaAssg <- aggregate(as.formula(paste(". ~", rankLevel, "+", groupLevel)), data=taxaAssgReads[,-1], FUN=sum)
  
  cat(matrixNames[expId], ":", rankLevel, "=", nrow(taxaAssg), ",", groupLevel, "=", length(unique(taxaAssg[,groupLevel])), ".\n")
  
  rownames(taxaAssg) <- taxaAssg[,1]
  taxaAssgPer <- prop.table(as.matrix(taxaAssg[,3:ncol(taxaAssg)]), 2)
  colSums(taxaAssgPer)
  
  data.dist=dist.JSD(taxaAssgPer)
  
  n.max=20
  nclusters=run.pam.cl(taxaAssgPer, data.dist, n.max)
  plotOptClusters(figPath, nclusters, n.max)
  
  # change k according to optimal-clusters-JSD.pdf
  k=9 #15, 12, 4
  data.cluster=pam.clustering(data.dist, k)
  pam.clustering(data.dist, k, as.vector=F)
  validateOptClusters(taxaAssgPer, data.dist, data.cluster)
  
  plotEnterotypes(file.path(workingPath, figDir), taxaAssgPer, data.dist, data.cluster, percent=0.01)
}
