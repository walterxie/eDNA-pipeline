# Reformatting of plot prioritisation figures

library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(reshape2)
library(corrgram)
library(NeatMap)
library(phyloseq)

theme_set(theme_bw(base_size=8))

setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/")
setwd("C:/Users/Andrew/PhD_folder/LBI_U8_OTUs_analyses/")
#setwd("J:/PhD/PhD Research/NZGL01401_analysis/LBI_U8_OTUs_analyses")

envdata.bysubplot <- read.table("LBI_environmental_data/LBI_all_env_data_by_subplot.txt", 
                                sep = "\t", header = TRUE, row.names = 1)
envdata.byplot <- read.table("LBI_environmental_data/LBI_all_env_data_by_plot.txt", 
                             sep = "\t", header = TRUE, row.names = 1)
envdata.byplot.sort <- envdata.byplot[order(envdata.byplot$Elevation),]
envdata.byplot.sort$Plot_elev <- paste0(rownames(envdata.byplot.sort), " (", round(envdata.byplot.sort$Elevation,0), "m)")

### Overall prioritisation data ###
d.overall <- read.table("Plot_prioritisation/Plot_priority_data_overall.txt", sep = "\t", header = TRUE)
d.phylo.overall <- read.table("Plot_prioritisation/LBI_phylo_alpha_plot_ranks_by_gene.txt", sep = "\t", header = TRUE)

d.phylo.overall$X <- NULL
d.phylo.overall$Elevation <- NULL
d.phylo.overall$ForestType <- NULL
d.phylo.overall$Variable <- "phylo"
d.phylo.overall$Mean <- rowMeans(d.phylo.overall[,2:7])
d.phylo.overall$Sd <-apply(d.phylo.overall[,2:7], 1, sd)
colnames(d.phylo.overall) <- gsub("\\.EUKARYOTA|\\.PROKARYOTA", "", colnames(d.phylo.overall))

d.overall <- rbind(d.overall, d.phylo.overall)
d.overall.l <- melt(d.overall[,1:8], id.vars = c("Variable", "Plot"), variable.name = "Gene", value.name = "Rank")

### Prioritisation data by group ###
d.group <- read.table("Plot_prioritisation/Plot_priority_data_by_group.txt", sep = "\t", header = TRUE)
d.phylo <- read.table("Plot_prioritisation/LBI_phylo_alpha_plot_ranks_by_group.txt", sep = "\t", header = TRUE)
d.phylo$X <- NULL
d.phylo$Elevation <- NULL
d.phylo$ForestType <- NULL
d.phylo$Variable <- "phylo"
d.phylo$Mean <- rowMeans(d.phylo[,3:15])
d.phylo$Sd <-apply(d.phylo[,3:15], 1, sd)
d.group <- rbind(d.group, d.phylo)
colnames(d.group) <- gsub("X", "", colnames(d.group))
colnames(d.group) <- gsub("FolCO1", "COI-650", colnames(d.group))
colnames(d.group) <- gsub("ShCO1", "COI-300", colnames(d.group))
colnames(d.group) <- gsub(".BACTERIA", " Bacteria", colnames(d.group))
colnames(d.group) <- gsub(".PROTISTS", " Protists", colnames(d.group))
colnames(d.group) <- gsub(".FUNGI", " Fungi", colnames(d.group))
colnames(d.group) <- gsub(".ANIMALIA", " Animals", colnames(d.group))

### Heatmap using NeatMap with column and row clustering ###
for(v in c("gamma0","gamma1","beta1","phylo")){
  a <- d.group[d.group$Variable == v, ]
  a$Elevation <- envdata.byplot$Elevation[match(a$Plot, rownames(envdata.byplot))]
  b <- a[, 3:16] # Get data columns
  rownames(b) <- paste0(a$Plot)#, " (", round(a$Elevation, 0), "m)")
  d <- dist(as.matrix(t(b))) # Euclidean distance
  h <- hclust(d)
  #c <- cor(b, method = "pearson")
  #h <- hclust(as.dist(1-c)) # Gives same clustering as Euclidean
  #plot(h)
  rownames(b) <- factor(rownames(b), levels = envdata.byplot.sort$Plot_elev, ordered = TRUE)
  p <- make.heatmap1(as.matrix(b), 
                column.method = "PCA", column.metric="euclidean", # Column ordering method
                row.method  ="complete", row.metric = "euclidean", 
                column.cluster.method = "complete", # Column clustering method
                row.cluster.method = "complete",
                column.labels = colnames(b), row.labels = rownames(b),
                column.label.size = 2.5, row.label.size = 2.5) +
                scale_fill_gradient2(low="#f46d43", mid="#ffffbf", high="#3288bd",
                                     midpoint = 14.5) + theme(legend.position="none") +
                xlab("") + ylab("") + #coord_flip() +
                theme(panel.grid = element_blank(), panel.border = element_blank(),
                      axis.text = element_blank(), axis.ticks = element_blank(),
                      plot.margin = unit(c(0,2,2,0),"cm"), strip.text.x = element_text(colour = "red"))
  
  
  # Turn off clipping
  gg1 <- ggplotGrob(p)
  gg1$layout[grepl("panel", gg1$layout$name), ]$clip <- "off"
  
  ### Get data for secondary plot ###
  gb = ggplot_build(p)
  #grid.draw(gb)
  
  #build$data[[4]]$y
  #build$data[[4]]$label
  p.ord <- as.data.frame(gb$data[[4]])
  #colnames(p.ord)[[1]] <- "Plot"
  p.ord <- p.ord[order(p.ord$y), ]
  p.ord$label <- factor(p.ord$label, levels = p.ord$label, ordered = TRUE)
  p.ord$Elevation <- envdata.byplot$Elevation[match(p.ord$label, rownames(envdata.byplot))]
  
  p2 <- ggplot(data = p.ord, aes(label, Elevation, group = 1)) + geom_line() + 
    coord_flip() + ylab("Elevation (m)") + xlab("") + 
    theme(axis.text.x = element_text(size = 8), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.margin = unit(c(0,0,0,-0.5), "cm"))
  gg2 <- ggplotGrob(p2)
  
  #build$data[[3]]$label <- gsub(" ", "\n", build$data[[3]]$label)
  #build$data[[3]]$angle <- "-45"
  #build$data[[3]]$vjust <- "1"
  #build$data[[3]]$hjust <- "1.5"
  
  #grid.draw(ggplot_gtable(build))

  maxHeight = grid::unit.pmax(gg1$heights[2:5], gg2$heights[2:5])
  gg1$heights[2:5] <- as.list(maxHeight)
  gg2$heights[2:5] <- as.list(maxHeight)
  
  # Don't know how to get clustered row order from heatmap (for adding elevation data series)
  #p2 <- ggplot(data = envdata.byplot, aes(Plot, Elevation, group = 1)) + geom_line() + coord_flip()
  pdf(file = paste0("LBI_", v, "_neatmap_by_group_clustered_x_and_y.pdf"), width = 8, height = 5, useDingbats = FALSE)
  plot(gt)
  grid.arrange(gg1, gg2, ncol=2, widths = c(10, 1))
  dev.off()
}  

### Heatmap using NeatMap with column clustering only ###
for(v in c("gamma0","gamma1","beta1","phylo")){
  a <- d.group[d.group$Variable == v, ]
  a$Elevation <- envdata.byplot$Elevation[match(a$Plot, rownames(envdata.byplot))]
  b <- a[, 3:16] # Get data columns
  rownames(b) <- paste0(a$Plot)#, " (", round(a$Elevation, 0), "m)")
  d <- dist(as.matrix(t(b))) # Euclidean distance
  h <- hclust(d)
  #c <- cor(b, method = "pearson")
  #h <- hclust(as.dist(1-c)) # Gives same clustering as Euclidean
  #plot(h)
  rownames(b) <- factor(rownames(b), levels = rownames(envdata.byplot.sort), ordered = TRUE)
  p <- make.heatmap1(as.matrix(b), 
                     column.method = "PCA", column.metric="euclidean", # Column ordering method
                     row.method  ="none", #row.metric = "euclidean", 
                     column.cluster.method = "complete", # Column clustering method
                     row.cluster.method = "none",
                     column.labels = colnames(b), row.labels = rownames(b),
                     column.label.size = 2.5, row.label.size = 2.5) +
    scale_fill_gradient2(low="#f46d43", mid="#ffffbf", high="#3288bd",
                         midpoint = 14.5) + theme(legend.position="none") +
    xlab("") + ylab("") + #coord_flip() +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          plot.margin = unit(c(0,2,2,0),"cm"))

  # Turn off clipping
  gt <- ggplotGrob(p)
  gt$layout[grepl("panel", gt$layout$name), ]$clip <- "off"
  pdf(file = paste0("LBI_", v, "_neatmap_by_group_clustered_x_only.pdf"), width = 8, height = 5, useDingbats = FALSE)
  plot(gt)
  dev.off()
} 
  
### Heatmaps using phyloseq ###
for(v in c("gamma0","gamma1","beta1","phylo")){
  a <- d.group[d.group$Variable == v, ]
  b <- a[, 3:16] # Get data columns
  rownames(b) <- a$Plot
  colnames(b) <- gsub(" ", "\n", colnames(b))
  
  phy <- phyloseq(otu_table(b, taxa_are_rows = FALSE))
  ph <- plot_heatmap(phy, method = "RDA", distance = "euclidean") 
  ph <- ph + scale_fill_gradient2(low="#f46d43", mid="#ffffbf", high="#3288bd", midpoint = 14.5) + 
       theme(legend.position="none", axis.text.x=element_text(size = 8, angle = -45, vjust = 1),
             axis.text.y = element_text(size = 8),
             panel.margin = unit(c(0,-0.5,0,0), "cm")) +
       ylab("Group") + xlab("Sample plot") + coord_flip()
  ph
  gg1 <- ggplotGrob(ph)
  
  ### Extract data for secondary plot ###
  build = ggplot_build(ph)
  p.ord <- as.data.frame(build$panel$ranges[[1]]$y.labels)
  colnames(p.ord)[[1]] <- "Plot"
  p.ord$Plot <- factor(p.ord$Plot, levels = p.ord$Plot, ordered = TRUE)
  p.ord$Elevation <- envdata.byplot$Elevation[match(p.ord$Plot, rownames(envdata.byplot))]
  
  p2 <- ggplot(data = p.ord, aes(Plot, Elevation, group = 1)) + geom_line() + 
        coord_flip() + ylab("Elevation (m)") + xlab("") + 
        theme(axis.text.x = element_text(size = 8), 
              axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              panel.margin = unit(c(0,0,0,-0.5), "cm"))
  gg2 <- ggplotGrob(p2)

  ### Make plots align ###
  #maxWidth = grid::unit.pmax(gg1$widths[2:5], gg2$widths[2:5])
  #gg1$widths[2:5] <- as.list(maxWidth)
  #gg2$widths[2:5] <- as.list(maxWidth)
  maxHeight = grid::unit.pmax(gg1$heights[2:5], gg2$heights[2:5])
  gg1$heights[2:5] <- as.list(maxHeight)
  gg2$heights[2:5] <- as.list(maxHeight)

  pdf(file = paste0("LBI_", v, "_phy-heatmap_by_group_clustered_x_and_y.pdf"), width = 8, height = 5, useDingbats = FALSE)
  grid.arrange(gg1, gg2, ncol=2, widths = c(10, 2))
  dev.off()

}
  
  
for(v in c("gamma0","gamma1","beta1","phylo")){
    a <- d.group[d.group$Variable == v, ]
    b <- a[, 3:16] # Get data columns
    rownames(b) <- a$Plot
    colnames(b) <- gsub(" ", "\n", colnames(b))
    
    phy <- phyloseq(otu_table(b, taxa_are_rows = FALSE))
    ph <- plot_heatmap(phy, method = "RDA", distance = "euclidean",
                       sample.order = rownames(envdata.byplot.sort)) 
    ph <- ph + scale_fill_gradient2(low="#f46d43", mid="#ffffbf", high="#3288bd", midpoint = 14.5) + 
      theme(legend.position="none", axis.text.x=element_text(size = 8, angle = -45, vjust = 1),
            axis.text.y = element_text(size = 8),
            panel.margin = unit(c(0,-0.5,0,0), "cm")) +
            ylab("Group") + xlab("Sample plot") + coord_flip()
    gg1 <- ggplotGrob(ph)
    
    ### Extract data for secondary plot ###
#     build = ggplot_build(ph)
#     p.ord <- as.data.frame(build$panel$ranges[[1]]$y.labels)
#     colnames(p.ord)[[1]] <- "Plot"
#     p.ord$Plot <- factor(p.ord$Plot, levels = p.ord$Plot, ordered = TRUE)
#     p.ord$Elevation <- envdata.byplot$Elevation[match(p.ord$Plot, rownames(envdata.byplot))]
#     
#     p2 <- ggplot(data = p.ord, aes(Plot, Elevation, group = 1)) + geom_line() + 
#       coord_flip() + ylab("Elevation (m)") + xlab("") + 
#       theme(axis.text.x = element_text(size = 8), 
#             axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#             panel.margin = unit(c(0,0,0,-0.5), "cm"))
#     gg2 <- ggplotGrob(p2)
    
    ### Make plots align ###
    #maxWidth = grid::unit.pmax(gg1$widths[2:5], gg2$widths[2:5])
    #gg1$widths[2:5] <- as.list(maxWidth)
    #gg2$widths[2:5] <- as.list(maxWidth)
#    maxHeight = grid::unit.pmax(gg1$heights[2:5], gg2$heights[2:5])
#    gg1$heights[2:5] <- as.list(maxHeight)
#    gg2$heights[2:5] <- as.list(maxHeight)
    
    pdf(file = paste0("LBI_", v, "_phy-heatmap_by_group_clustered_x_only.pdf"), width = 8, height = 5, useDingbats = FALSE)
#    grid.arrange(gg1, gg2, ncol=2, widths = c(10, 2))
    plot(ph)
    dev.off()
  
}
  
  
##########################################################  
### Old stuff ###  
  
  
                
  m <- melt(a[, !(grepl("Mean|Sd", colnames(a)))], id.vars = c("Variable","Plot"), variable.name = "group", value.name = "rank")
  m$Plot <- factor(m$Plot, levels = rownames(envdata.byplot.sort), ordered = TRUE)
  m$group <- factor(m$group, levels = h$labels[c(h$order)], ordered = TRUE) 
  
  #d.group.l <- d.group.l[!(d.group.l$Gene == "COI-650" & d.group.l$Taxon == "Fungi"), ] # Drop COI-650 fungi data
  p <- ggplot(m, aes(x = group, y = Plot)) + geom_tile(aes(fill = rank)) +
       scale_fill_gradient2(low="#f46d43", mid="#ffffbf", high="#3288bd",
                            midpoint = 14.5) + theme(legend.position="none") +
       xlab("") + ylab("") +
       #facet_grid(.~Taxon, scales = "free", space = "free", drop=TRUE) + 
       theme(strip.background = element_blank())
    # Saved at 8 x 4 in (landscape)
   ggsave(p, file = paste0("LBI_", v, "_plot_ranks_by_group_recol_v2.pdf"), height = 4, width = 8, useDingbats = FALSE)

}




d.group.l <- melt(d.group[,1:16], id.vars = c("Variable", "Plot"), variable.name = "Gene", value.name = "Rank")

### Heatmap plots ###
d.overall.l$Elevation <- envdata.byplot$Elevation[match(d.overall.l$Plot, rownames(envdata.byplot))]
d.group.l$Elevation <- envdata.byplot$Elevation[match(d.group.l$Plot, rownames(envdata.byplot))]

d.overall.l$Gene <- gsub("X", "", d.overall.l$Gene)
d.overall.l$Gene <- gsub("FolCO1", "COI-650", d.overall.l$Gene)
d.overall.l$Gene <- gsub("ShCO1", "COI-300", d.overall.l$Gene)
d.overall.l$Gene <- factor(d.overall.l$Gene, levels = c("16S", "18S", "26S", "ITS", "COI-300", "COI-650", "Mean"), ordered = TRUE)
d.overall.l$Plot <- factor(d.overall.l$Plot, levels = rownames(envdata.byplot.sort), ordered = TRUE)

d.group.l$Gene <- gsub("X", "", d.group.l$Gene)
d.group.l$Gene <- gsub("FolCO1", "COI-650", d.group.l$Gene)
d.group.l$Gene <- gsub("ShCO1", "COI-300", d.group.l$Gene)
d.group.l$Gene <- gsub(".BACTERIA", " Bacteria", d.group.l$Gene)
d.group.l$Gene <- gsub(".PROTISTS", " Protists", d.group.l$Gene)
d.group.l$Gene <- gsub(".FUNGI", " Fungi", d.group.l$Gene)
d.group.l$Gene <- gsub(".ANIMALIA", " Animals", d.group.l$Gene)
d.group.l$Taxon <- sapply(strsplit(as.character(d.group.l$Gene), "\\s"), "[[", 2)
#d.group.l$Gene <- sapply(strsplit(as.character(d.group.l$Gene), "\\s"), "[[", 1)
d.group.l$Gene <- gsub(" Bacteria| Protists| Fungi| Animals", "", d.group.l$Gene)
d.group.l$Gene <- factor(d.group.l$Gene, levels = c("16S","18S","26S","ITS","COI-300","COI-650"), ordered = TRUE)
# d.group.l$Gene <- gsub(" ", "\n", d.group.l$Gene) # For line breaks in axis labels
# d.group.l$Gene <- factor(d.group.l$Gene, levels = c("16S\nBacteria","18S\nProtists","26S\nProtists","COI-300\nProtists", 
#                                                     "COI-650\nProtists","18S\nFungi","26S\nFungi","ITS\nFungi","COI-300\nFungi",
#                                                     "COI-650\nFungi","18S\nAnimals","26S\nAnimals","COI-300\nAnimals","COI-650\nAnimals",
#                                                     "Mean"), ordered = TRUE)
d.group.l$Plot <- factor(d.group.l$Plot, levels = rownames(envdata.byplot.sort), ordered = TRUE)
d.group.l$Taxon <- factor(d.group.l$Taxon, levels = c("Bacteria", "Protists", "Fungi", "Animals"), ordered = TRUE)

for(v in unique(d.overall.l$Variable)){
#d.overall.l.2 <- d.overall.l[d.overall.l$Gene != "Mean",]
  p <- ggplot(d.overall.l[d.overall.l$Variable == "gamma0", ], 
       aes(x = Gene, y = Plot)) + geom_tile(aes(fill = Rank)) +
      #scale_fill_gradientn(colors = c("White", "SteelBlue"), breaks = c(1, 7, 14, 21, 28))
      scale_fill_gradient2(low="#f46d43", mid="#ffffbf", high="#3288bd",
                           midpoint = 14.5) + theme(legend.position="none") # Exported at 4 x 5 inches
  ggsave(p, file = paste0("LBI_", v, "_plot_ranks_by_gene_recol.pdf"), height = 4, width = 8, useDingbats = FALSE)
}  

d.group.l <- d.group.l[!(d.group.l$Gene == "COI-650" & d.group.l$Taxon == "Fungi"), ] # Drop COI-650 fungi data
for(v in unique(d.group.l$Variable)){
  p <- ggplot(d.group.l[d.group.l$Variable == v, ], aes(x = Gene, y = Plot)) + geom_tile(aes(fill = Rank)) +
    #scale_fill_gradientn(colors = c("White", "SteelBlue"), breaks = c(1, 7, 14, 21, 28))
    scale_fill_gradient2(low="#f46d43", mid="#ffffbf", high="#3288bd",
                         midpoint = 14.5) + theme(legend.position="none") +
    facet_grid(.~Taxon, scales = "free", space = "free", drop=TRUE) + 
    theme(strip.background = element_blank())
    # Saved at 8 x 4 in (landscape)
  ggsave(p, file = paste0("LBI_", v, "_plot_ranks_by_group_recol_v2.pdf"), height = 4, width = 8, useDingbats = FALSE)
}


### Plot prioritisation correlations ###
d.overall.x <- dcast(d.overall.l, Plot ~ Gene + Variable, value.var = c("Rank"))
colnames(d.overall.x) <- gsub("_", "\n", colnames(d.overall.x))
pdf(file = paste0("Overall_gene_prioritisation_correlations_2.pdf"), height = 15, width = 15)
plot(d.overall.x[ ,!(grepl("Plot", colnames(d.overall.x)))], gap = 0, 
     lower.panel = panel.smooth, upper.panel = panel.conf, 
     cex.axis = 0.75, cex.cor = 0.9, col.smooth = "purple") 
dev.off()

vars <- c("gamma0","gamma1","beta1","phylo")
for(v in vars){
  d.overall.x <- dcast(d.overall.l[grepl(v, d.overall.l$Variable), ], Plot ~ Gene, value.var = c("Rank"))
  colnames(d.overall.x) <- gsub("_", "\n", colnames(d.overall.x))
  pdf(file = paste0(v, "_gene_prioritisation_correlations_2.pdf"), height = 4, width = 4)
  plot(d.overall.x[ ,!(grepl("Plot", colnames(d.overall.x)))], gap = 0, 
       lower.panel = panel.smooth, upper.panel = panel.conf, 
       cex.axis = 0.75, cex.cor = 0.9, col.smooth = "purple") 
  dev.off()
}

vars <- c("16S","18S","26S","ITS","COI-300","COI-650")
for(v in vars){
  d.overall.x <- dcast(d.overall.l[grepl(v, d.overall.l$Gene), ], Plot ~ Variable, value.var = c("Rank"))
  colnames(d.overall.x) <- gsub("_", "\n", colnames(d.overall.x))
  pdf(file = paste0(v, "_var_prioritisation_correlations_2.pdf"), height = 4, width = 4)
  plot(d.overall.x[ ,!(grepl("Plot", colnames(d.overall.x)))], gap = 0, 
       lower.panel = panel.smooth, upper.panel = panel.conf, 
       cex.axis = 0.75, cex.cor = 0.9, col.smooth = "purple") 
  dev.off()
}

vars <- c("gamma0","gamma1","beta1","phylo")
for(v in vars){
  d.group.x <- dcast(d.group.l[grepl(v, d.group.l$Variable), ], Plot ~ Gene + Taxon, value.var = c("Rank"))
  colnames(d.group.x) <- gsub(paste0(v, "_"), "", colnames(d.group.x))
  colnames(d.group.x) <- gsub("_", "\n", colnames(d.group.x))
  pdf(file = paste0(v, "_group_prioritisation_correlations_2.pdf"), height = 10, width = 10)
  plot(d.group.x[ ,!(grepl("Plot", colnames(d.group.x)))], gap = 0, 
       lower.panel = panel.smooth, upper.panel = panel.conf, 
       cex.axis = 0.75, cex.cor = 0.9, col.smooth = "purple")  
  dev.off()
}

vars <- c("16S","18S","26S","ITS","COI-300","COI-650")
for(v in vars){
  d.group.x <- dcast(d.group.l[grepl(v, d.group.l$Gene), ], Plot ~ Variable + Taxon, value.var = c("Rank"))
  colnames(d.group.x) <- gsub(paste0(v, "_"), "", colnames(d.group.x))
  colnames(d.group.x) <- gsub("_", "\n", colnames(d.group.x))
  z <- ncol(d.group.x)-1
  pdf(file = paste0(v, "_group_prioritisation_correlations_2.pdf"), height = z/1.5, width = z/1.5)
  plot(d.group.x[ ,!(grepl("Plot", colnames(d.group.x)))], gap = 0, 
       lower.panel = panel.smooth, upper.panel = panel.conf, 
       cex.axis = 0.75, cex.cor = 0.9, col.smooth = "purple")  
  dev.off()
}

