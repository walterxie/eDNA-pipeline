library(RColorBrewer)

if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("tableFile")) stop("table file is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("taxa.group")) taxa.group="assigned"
if(!exists("isPlot")) isPlot = FALSE # by subplot

#source("Modules/init.r")
source("Modules/RarefactionTable.R")
source("Modules/RarefactionPlot.R")

n <- length(matrixNames) 

cat("\nAnalysis: plot rarefaction table per sample. \n")

env <- getSampleMetaData(isPlot)
env[,"ForestType"] <- gsub(":.*", "", env[,"ForestType"], ignore.case = T)
env[,"ForestType"] <- gsub("x", "unknown", env[,"ForestType"], ignore.case = T)

for (expId in 1:(n-1)) {
  ### alpha0 ###
  div="alpha0"
  rare.df <- getRarefactionTableTaxa(expId, isPlot, rmSingleton, taxa.group, div)
  
  if (is.null(rare.df)) {
    cat("\nSkip", taxa.group, "subset from", matrixNames[expId], "because of no", div, "rarefaction table file.\n")
    next
  }
  
  fname <- paste("rf", div, matrixNames[expId], postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-")
  fileStem <- file.path(workingPath, figDir, fname)
  plotRarefactionTable(rare.df, env, fileStem, ylab=div)
  
  ### alpha1 ###
  div="alpha1"
  rare.df <- getRarefactionTableTaxa(expId, isPlot, rmSingleton, taxa.group, div)
  
  if (is.null(rare.df)) {
    cat("\nSkip", taxa.group, "subset from", matrixNames[expId], "because of no", div, "rarefaction table file.\n")
    next
  }
  
  fname <- paste("rf", div, matrixNames[expId], postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-")
  fileStem <- file.path(workingPath, figDir, fname)
  plotRarefactionTable(rare.df, env, fileStem, ylab=div)
  
}
