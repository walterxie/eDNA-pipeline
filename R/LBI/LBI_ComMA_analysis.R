library(ComMA)
library(picante)
theme_set(theme_bw(base_size=8))

#setwd("C:/Documents and Settings/Andrew/Desktop/LBI_miseq_analysis_stuff/")
setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/ComMA_stuff/")

### Load environmental data ###
envdata.bysubplot <- read.table("LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_subplot.txt", 
                                sep = "\t", header = TRUE, row.names = 1)
envdata.byplot <- read.table("LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_plot.txt", 
                             sep = "\t", header = TRUE, row.names = 1)



### https://github.com/walterxie/eDNA-pipeline/blob/hauturu/R/allPlotPriorPDDiversity.r ###
if(!exists("tableFile")) stop("table file is missing !")
if(!exists("figDir")) stop("figure folder name is missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("otuThr")) otuThr = 97


######## heatmap #######

matrixNames <-  c("16S", "18S", "26S", "ITS", "ShCO1", "FolCO1")
#taxaGroups <- c("BACTERIA", "FUNGI", "PROTISTS", "ANIMALIA") # default rank=kingdom
taxaGroups <- c("PROKARYOTA","EUKARYOTA")
env <- getSampleMetaData(TRUE)

for (taxag in taxaGroups) {
  for (matrn in matrixNames) {
#    if (matrn != "16S" && taxag == "BACTERIA")
#      next
     if (matrn != "16S" && taxag == "PROKARYOTA") 
       next
     if (matrn == "16S" && taxag == "EUKARYOTA")
       next
    
    t.communityMatrix <- getCommunityMatrixT(matrix.name=matrn, isPlot=TRUE, minAbund=2, minRich=200, taxa.group="all")
    
    if (!is.null(t.communityMatrix)) {
      phylo.tree <- getPhyloTree(matrn, taxag)
      if (!is.null(phylo.tree)) {
        combined <- match.phylo.comm(phylo.tree, t.communityMatrix)
        pd.alpha <- phylo.alpha(combined$comm, combined$phy)
        #pd.alpha <- phylo.alpha(t.communityMatrix, phylo.tree)
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
colnames(rank.elv.PD)[1] <- "Plot"

write.table(rank.elv.PD, file = "LBI_phylo_alpha_plot_ranks_by_group.txt", sep = "\t", quote = F, col.names = NA)
write.table(rank.elv.PD, file = "LBI_phylo_alpha_plot_ranks_by_gene.txt", sep = "\t", quote = F, col.names = NA)

### Make a heatmap (not using ComMA method) ###
m <- melt(rank.elv.PD, id = c("Plot","Elevation","ForestType"), value.name = "rank")
m$Gene <- gsub("\\.\\w*", "", m$variable)
m$Gene <- gsub("FolCO1","COI-650", m$Gene)
m$Gene <- gsub("ShCO1","COI-300", m$Gene)
m$Gene <- factor(m$Gene, levels = c("16S","18S","26S","ITS","COI-300","COI-650"), ordered = TRUE)
m <- m[order(m$Elevation),]
m$Plot <- factor(m$Plot, levels = m$Plot, ordered = TRUE)
ggplot(m, aes(x = gene, y = Plot, fill = rank)) + geom_tile() +
  scale_fill_gradient2(low="#f46d43", mid="#ffffbf", high="#3288bd",
                       midpoint = 14.5) + theme(legend.position="none") # Exported at 4 x 5 inches


m <- melt(rank.elv.PD, id = c("Plot","Elevation","ForestType"), value.name = "rank")
m <- m[!(m$variable == "FolCO1.FUNGI"), ] # Drop COI-650 fungi data
m$Gene <- sapply(strsplit(as.character(m$variable), "\\."), "[[", 1)
m$Gene <- gsub("FolCO1","COI-650", m$Gene)
m$Gene <- gsub("ShCO1","COI-300", m$Gene)
m$Gene <- factor(m$Gene, levels = c("16S","18S","26S","ITS","COI-300","COI-650"), ordered = TRUE)
m$Taxon <- sapply(strsplit(as.character(m$variable), "\\."), "[[", 2)
m$Taxon <- gsub("BACTERIA", "Bacteria", m$Taxon)
m$Taxon <- gsub("PROTISTS", "Protists", m$Taxon)
m$Taxon <- gsub("FUNGI", "Fungi", m$Taxon)
m$Taxon <- gsub("ANIMALIA", "Animals", m$Taxon)
m$Taxon <- factor(m$Taxon, levels = c("Bacteria", "Protists", "Fungi", "Animals"), ordered = TRUE)
m$Plot <- factor(m$Plot, levels = m$Plot, ordered = TRUE)
p <- ggplot(m, aes(x = Gene, y = Plot, fill = rank)) + geom_tile() +
  scale_fill_gradient2(low="#f46d43", mid="#ffffbf", high="#3288bd",
                       midpoint = 14.5) + theme(legend.position="none") +
  facet_grid(.~Taxon, scales = "free", space = "free", drop=TRUE) + 
  theme(strip.background = element_blank()) # Saved at 8 x 4 in (landscape)

ggsave(p, file = paste0("LBI_phylo_alpha_plot_ranks_by_group_recol_v2.pdf"), height = 4, width = 8, useDingbats = FALSE)

### ComMA heatmap ###
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