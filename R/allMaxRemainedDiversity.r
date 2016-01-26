
library(gplots)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(xtable)


if(!exists("tableFile")) stop("table file is missing !")
if(!exists("figDir")) stop("figure folder name is missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("otuThr")) otuThr = 97

#source("Modules/init.R", local=TRUE)
source("Modules/MaximizeDiversity.R")

env <- getSampleMetaData(TRUE)
env[,"ForestType"] <- gsub(":.*", "", env[,"ForestType"], ignore.case = T)
env[,"ForestType"] <- gsub("x", "unknown", env[,"ForestType"], ignore.case = T)

######## heatmap #######
divs <- c("gamma-1", "beta-1")
subTitles <- c(expression(paste("Ranked by "^"1", D[gamma])),expression(paste("Ranked by "^"1", D[beta])))

taxa.group="assigned"
for (i in 1:length(divs)) {
  # ranks, last 2 col are mean and sd
  allRanks <- getMaxRemainedDiversityRank(divs[i], taxa.group) 
  rank.elv <- merge(allRanks, env[,c("Elevation","ForestType")], by = "row.names")
  rank.elv <- rank.elv[order(rank.elv[,"Elevation"]),]
  colnames(rank.elv)[1] <- "Samples"

  fname <- paste("max-remained", divs[i], taxa.group,"heatmap.pdf", sep = "-")
  plotHeatmapMaxRemDivRank(rank.elv, fname, title=subTitles[i])
  
  rownames(rank.elv) <- rank.elv[,"Samples"]
  rank.elv <- rank.elv[,-c(1,ncol(rank.elv)-1,ncol(rank.elv))]
  printXTable(rank.elv, caption = paste("Ranking sampling plots by removing plots sequentially so as to 
			 minimize the loss of overall", divs[i], "diversity among the remaining plots. 
	     1 is the most important and removed at the last, 10 is the least important and removed in the beginning."), 
       label = paste0("tab:maxRemained", divs[i]), file=tableFile)
}

taxaGroups <- c("ANIMALIA", "FUNGI", "PROTISTS")
for (i in 1:length(divs)) {
  # hard code
  allRanks <- getMaxRemainedDiversityRank(divs[i], taxa.group="BACTERIA") 
  allRanksGourp <- data.frame(row.names = rownames(allRanks), stringsAsFactors = F)
  allRanksGourp[,"16S.BACTERIA"] <- allRanks[,"16S"]
  
  for (taxag in taxaGroups) {
      taxa.group <- taxag
      allRanks <- getMaxRemainedDiversityRank(divs[i], taxa.group=taxa.group) 
      colnames(allRanks) <- paste(colnames(allRanks), taxa.group, sep = ".")
      
      allRanksGourp <- merge(allRanksGourp, allRanks[,-c(ncol(allRanks)-1,ncol(allRanks))], by = "row.names")
      rownames(allRanksGourp) <- allRanksGourp[,"Row.names"]
      allRanksGourp <- allRanksGourp[,-1]
  }

  rank.elv <- merge(allRanksGourp, env[,c("Elevation","ForestType")], by = "row.names")
  rank.elv <- rank.elv[order(rank.elv[,"Elevation"]),]
  colnames(rank.elv)[1] <- "Samples"
  
  fname <- paste("max-remained", divs[i], "group-heatmap.pdf", sep = "-")
  plotHeatmapMaxRemDivRank(rank.elv, fname, title=subTitles[i])
  
  rownames(rank.elv) <- rank.elv[,"Samples"]
  rank.elv <- rank.elv[,-c(1,ncol(rank.elv)-1,ncol(rank.elv))]
  printXTable(rank.elv, caption = paste("Ranking sampling plots by removing plots sequentially so as to 
              minimize the loss of overall", divs[i], "diversity among the remaining plots. 
              1 is the most important and removed at the last, 10 is the least important and removed in the beginning."), 
              label = paste0("tab:maxRemained", divs[i], "Group"), file=tableFile)
}

