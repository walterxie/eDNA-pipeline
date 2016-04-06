
library(ComMA)


if(!exists("tableFile")) stop("table file is missing !")
if(!exists("figDir")) stop("figure folder name is missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("otuThr")) otuThr = 97


######## heatmap #######

matrixNames <-  c("16S", "18S", "26S", "ITS", "ShCO1", "FolCO1")
taxaGroups <- c("BACTERIA", "FUNGI", "PROTISTS", "ANIMALIA") # default rank=kingdom
env <- getSampleMetaData(TRUE)

for (taxag in taxaGroups) {
  for (matrn in matrixNames) {
    if (matrn != "16S" && taxag == "BACTERIA") 
      next
    
    t.communityMatrix <- getCommunityMatrixT(matrn, isPlot=TRUE, minAbund=2, minRich=600, taxa.group=taxag)
    
    if (!is.null(t.communityMatrix)) {
      phylo.tree <- getPhyloTree(matrn, taxag)
      if (!is.null(phylo.tree)) {
        pd.alpha <- phylo.alpha(t.communityMatrix, phylo.tree)
        pd.alpha <- pd.alpha[order(pd.alpha$PD, decreasing = T), ]
        pd.alpha$rank.PD <- 1:nrow(pd.alpha)
        pd.alpha <- pd.alpha[order(pd.alpha$SR, decreasing = T), ]
        pd.alpha$rank.SR <- 1:nrow(pd.alpha)
        pd.alpha <- pd.alpha[order(rownames(pd.alpha)), ]
        
        if (matrn == "16S") {
          rank.PD.df <- data.frame(row.names = rownames(pd.alpha))
          rank.SR.df <- data.frame(row.names = rownames(pd.alpha))
        } 
        colnam <- paste(matrn, taxag, sep = ".")
        rank.PD.df[,colnam] <- pd.alpha$rank.PD
        rank.SR.df[,colnam] <- pd.alpha$rank.SR
      }
    }
  }
}

rank.elv.PD <- merge(rank.PD.df, env[,c("Elevation","ForestType")], by = "row.names")
rank.elv.PD <- rank.elv.PD[order(rank.elv.PD[,"Elevation"]),]
colnames(rank.elv.PD)[1] <- "plot"

gg.plot <- ggHeatmap(df=rank.elv.PD[,-c(ncol(rank.elv.PD)-1,ncol(rank.elv.PD))], 
                     melt.id="plot", title="By Phylogenetic Alpha Diversity")
pdfGgplot(gg.plot, fig.path="plot-prior-pd-heatmap.pdf")

printXTable(rank.elv.PD[,-c(ncol(rank.elv.PD)-1,ncol(rank.elv.PD))], 
            caption = paste("Ranking sampling plots by phylogenetic alpha diversity proposed by Faith (1992). 
              1 is the most important and its plot has the biggest phylogenetic diversity."), 
            label = paste0("tab:plotPriorPD"), file=tableFile)

rank.elv.SR <- merge(rank.SR.df, env[,c("Elevation","ForestType")], by = "row.names")
rank.elv.SR <- rank.elv.SR[order(rank.elv.SR[,"Elevation"]),]
colnames(rank.elv.SR)[1] <- "plot"

gg.plot <- ggHeatmap(df=rank.elv.SR[,-c(ncol(rank.elv.SR)-1,ncol(rank.elv.SR))], 
                     melt.id="plot", title="By OTU Richness")
pdfGgplot(gg.plot, fig.path="plot-prior-sr-heatmap.pdf")

printXTable(rank.elv.SR[,-c(ncol(rank.elv.SR)-1,ncol(rank.elv.SR))], 
            caption = paste("Ranking sampling plots by OTU richness. 
              1 is the most important and its plot has the biggest OTU richness."), 
            label = paste0("tab:plotPriorSR"), file=tableFile)

