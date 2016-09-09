
library(igraph)
library(GGally)
library(reshape2)
library(ggplot2)
library(data.table)
theme_set(theme_bw(base_size=8))
library(RColorBrewer)
colors <- brewer.pal(4, "Spectral")
#colors <- brewer.pal(4, "Paired")
#colors <- brewer.pal(4, "Set3")
pal <- colorRampPalette(colors) 
source("H:/My Documents/PhD Research PFR folder/R stuff/R_LBI_Miseq/IWantHue.R")

#for(l in c("a","b","c","d")){
for(l in c("16S-18S26SCOI","16S-18SITSCOI")){
# Path to network data
#setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/Network_analysis/Networks_by_OTUs_various_BLAST_thresholds_percent_all_genes_v3/")
setwd(paste0("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/Network_analysis/Networks_by_OTUs_all_genes_veg_v3_phy_20-80_16S-18S26SCOI/"))#, l, "/"))
#setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/Network_analysis/Networks_by_OTUs_various_BLAST_thresholds_percent_all_genes_v4/")
#setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/Network_analysis/Networks_by_OTUs_various_BLAST_thresholds_percent_all_genes_veg_v4/")

e.file <- Sys.glob("*edges_data_90_v3_20-80.txt")
v.file <- Sys.glob("*vertices_data_90_v3_20-80.txt")

label <- gsub("vertices_data_90_v3_20-80.txt", "", v.file)
outfile <- paste0(label, "communities_summary_90_v3_20-80.txt")

### Recreate network structure ###
verts.mb <- read.table(v.file, sep = "\t", header = TRUE)
edges.mb <- read.table(e.file, sep = "\t", header = TRUE, row.names = 1)
g <- graph.data.frame(edges.mb, directed=FALSE, vertices=verts.mb)  

### Extract/recreate vegetation network structure only? ###
#verts.mb.2 <- verts.mb[verts.mb$Gene == "Vegetation", ]
#edges.mb.2 <- edges.mb[edges.mb$from.k == "Vegetation" & edges.mb$to.k == "Vegetation", ]
#g <- graph.data.frame(edges.mb.2, directed=FALSE, vertices=verts.mb.2)    

### Community detection using fast greedy algorithm ###
E(g)$weight <- abs(E(g)$weight) # Community detection doesn't work with negative values
grr <- fastgreedy.community(g)
length(grr)
sizes(grr)
modularity(grr)
#communities(grr)

#cat(paste("\n\n", v.file," fastgreedy communities\n"), file = outfile, append = T)
#write.table(as.data.frame(sizes(grr)), quote = FALSE, col.names = NA, file = outfile, append = TRUE)

# grr <- walktrap.community(g, steps = 5)
# length(grr)
# sizes(grr)
# modularity(grr)
# communities(grr)

# grr <- spinglass.community(g) #(slow)
# length(grr)
# sizes(grr)
# modularity(grr)
# communities(grr)

# grr <- infomap.community(g)
# length(grr)
# sizes(grr)
# modularity(grr)
# communities(grr)


### Plot each sub-network ###
# palette.gene = c("16S" = "#e46d1a", "18S" = "#377eb8",
#                  "26S" = "#4daf4a", "ITS" = "#984ea3", "COI-300" = "#e41a1c", 
#                  "COI-650" = "#f781bf", "Vegetation" = "#ffff33")
palette.gene = c("16S" = "#e46d1a", "18S" = "#377eb8",
                 "26S" = "#4daf4a", "ITS" = "#4daf4a", "COI-300" = "#f781bf", 
                 "COI-650" = "#e41a1c", "Vegetation" = "#ffff33")

#pdf(file=paste0(label, "fg_comm_all_graphs.pdf"), width=15/2.5, height=15/2.5, useDingbats = FALSE)
#plist <- list()
n <- 1
all.verts <- data.frame()
for(com in communities(grr)){
  #str(com)
  verts.com <- verts.mb[verts.mb$X %in% com, ]
  verts.com$comm <- n
  # Subset matching edges
  selectedRows <- (edges.mb$from %in% verts.com$X & edges.mb$to %in% verts.com$X)
  edges.com <- edges.mb[selectedRows,]
  g.com <- graph.data.frame(edges.com, directed=FALSE, vertices=verts.com)
  summary(g.com)
  
  # Vertex statistics
  V(g.com)$degree <- igraph::degree(g.com)
  g.com.uw <- g.com
  E(g.com.uw)$weight <- abs(E(g.com.uw)$weight) 
  g.com$vertices$betw <- igraph::betweenness(g.com.uw)
  betw <- igraph::betweenness(g.com.uw, normalized = TRUE)
  betw.norm <- (betw - min(betw))/(max(betw) - min(betw))
  g.com$vertices$betw_norm <- betw.norm
  # Edge statistics
  g.com$edges$betw <- igraph::edge.betweenness(g.com.uw) # Using weighted graph gives error
  #summary(g.com)
  
#   net <- as_data_frame(g.com, what = "both")
#   write.table(net$vertices, file=paste0(label, "_fg_comm_", n, "_vertices_data.txt"), 
#               sep = "\t", quote = FALSE, col.names = NA)
#   write.table(net$edges, file=paste0(label, "_fg_comm_", n, "_edges_data.txt"), 
#               sep = "\t", quote = FALSE, col.names = NA)
  
  # Print community graphs
  lty <- sapply(E(g.com)$weight, function(x) ifelse(x < 0, 2, 1))
#   g1 <- ggnet2(g.com, node.color = V(g.com)$Gene, palette = palette.gene, node.alpha = 0.5, 
#          node.label = V(g.com)$Phylum, node.size = V(g.com)$degree, 
#          label.size = 2, label.alpha = 0.5,
#          edge.size = abs(E(g.com)$weight*5), edge.alpha = 0.5, edge.lty = lty,
#          #edge.color = sign.mb,
#          edge.color = c("color", "grey50"),
#          legend.position = "none")
#   g2 <- ggnet2(g.com, node.color = V(g.com)$Gene, palette = palette.gene, node.alpha = 0.5, 
#          node.label = V(g.com)$Class, node.size = V(g.com)$degree, 
#          label.size = 2, label.alpha = 0.5,
#          edge.size = abs(E(g.com)$weight*5), edge.alpha = 0.5, edge.lty = lty,
#          #edge.color = sign.mb,
#          edge.color = c("color", "grey50"),
#          legend.position = "none")
#   g3 <- ggnet2(g.com, node.color = V(g.com)$Gene, palette = palette.gene, node.alpha = 0.5, 
#          node.label = V(g.com)$Genus, node.size = V(g.com)$degree, 
#          label.size = 2, label.alpha = 0.5,
#          edge.size = abs(E(g.com)$weight*5), edge.alpha = 0.5, edge.lty = lty,
#          #edge.color = sign.mb,
#          edge.color = c("color", "grey50"),
#          legend.position = "none")
  g4 <- ggnet2(g.com, node.color = V(g.com)$Gene, palette = palette.gene, node.alpha = 0.5, 
         node.size = V(g.com)$degree, 
         label.size = 2, label.alpha = 0.5,
         edge.size = abs(E(g.com)$weight*5), edge.alpha = 0.5, edge.lty = lty,
         #edge.color = sign.mb,
         edge.color = c("color", "grey50"),
         legend.position = "none")
 
  ggsave(paste0(label, "_fg_comm_", n, ".jpeg"), plot = g4, width = 15, height = 15, units = c("cm"), dpi = 200) 
   
  #pdf(file=paste0(label, "_fg_comm_", n, "_graphs.pdf"), width=15/2.5, height=15/2.5, useDingbats = FALSE)
  #plot(g1)
  #plot(g2)
  #plot(g3)
#  plot(g4)
  #dev.off()

#   # Output summary data
#   cat(paste("\nCommunity no.", n, "summary statistics\n"), file = outfile, append = TRUE)
#   cat(paste("\nVertices:", vcount(g.com), "Edges:", ecount(g.com), "\n"), file = outfile, append = TRUE)
#   cat(paste("\nGraph density:", graph.density(g.com)), file = outfile, append = TRUE)
#   cat(paste("\nTransitivity:", transitivity(g.com, type="global")), file = outfile, append = TRUE)
#   cat("\nDegree distribution:", file = outfile, append = TRUE)
#   write.table(as.data.frame(degree.distribution(g.com)), sep = "\t", quote = FALSE, col.names = NA, file = outfile, append = TRUE)
#   
#   dt <- data.table(verts.com)
#   zz <- dt[, .N, by=.(Kingdom, Gene, Phylum, Class, Order)]
#   cat("\nComposition summary:", file = outfile, append = TRUE)
#   write.table(as.data.frame(zz), sep = "\t", quote = FALSE, col.names = NA, file = outfile, append = TRUE)
# 
   all.verts <- rbind(all.verts, verts.com) 
  
  n <- n + 1
  
}

#dev.off()

### Output table with all community data ###
write.table(all.verts, file=paste0(label, "fg_all_communities_vertices_data.txt"), 
            sep = "\t", quote = FALSE, col.names = NA)

### Make plot of kingdom-level community composition ###
all.verts$Gene <- factor(all.verts$Gene, levels = c("16S","18S","26S","ITS","COI-300","COI-650","Vegetation"), ordered = TRUE)
library(ggplot2)
p1 <- ggplot(all.verts, aes(x = Gene, fill = Kingdom)) + geom_bar(position = "stack") + 
      scale_fill_brewer(palette = "Set2") +
      #scale_fill_manual(values = pal(7)) +  
      xlab("Amplicon or vegetation dataset") + 
      ylab(paste0("Number of OTUs")) +
      theme(strip.background = element_blank(), plot.title = element_text(size = 9),
        panel.grid = element_blank(), axis.text.x=element_text(angle = -90, hjust = 0)) +
      facet_wrap(~ comm)
ggsave(paste0(label, "fg_comm_composition_barplot.pdf"), 
       plot = p1, height = 15, width = 18, units = c("cm"), useDingbats = FALSE)
  

### Make a plot of phylum-level community composition ###
tax_ref <- read.table("H:/My Documents/PhD Research PFR folder/Database stuff/New_taxonomy_from_PLOSONE_2015_fixed.txt", 
                      header = TRUE, sep = "\t", quote = "", comment.char = "")
tax_ref$Phylum <- gsub("(\\s\\[=.*\\])", "", tax_ref$Phylum)

d <- read.table("LBI_all_genes_veg_ggnet_by_OTUs_fg_all_communities_vertices_data.txt", sep = "\t", header = T)

d$pos <- match(d$Phylum, tax_ref$Phylum)
d$Phylum <- with(d, reorder(Phylum, pos))
d$comm <- factor(d$comm)
d$Gene <- factor(d$Gene, levels = c("16S","18S","26S","ITS","COI-300","COI-650","Vegetation"), ordered = TRUE)

# ggplot(d[d$Gene == "16S",], aes(x = comm, fill = Phylum)) + geom_bar() 
# ggplot(d[d$Gene == "18S",], aes(x = comm, fill = Phylum)) + geom_bar() 
# ggplot(d[d$Gene == "26S",], aes(x = comm, fill = Phylum)) + geom_bar() 
# ggplot(d[d$Gene == "ITS",], aes(x = comm, fill = Phylum)) + geom_bar() 
# ggplot(d[d$Gene == "COI-300",], aes(x = comm, fill = Phylum)) + geom_bar() 
# ggplot(d[d$Gene == "COI-650",], aes(x = comm, fill = Phylum)) + geom_bar() 
# 
# ggplot(d[grepl("6|12|16", d$comm) & d$Gene == "16S",], aes(x = comm, fill = Phylum)) + geom_bar(position = "fill") + facet_grid(.~ Gene)

ppal <- iwanthue(length(unique(d$Phylum)))

p2 <- ggplot(d[d$Gene != "Vegetation",], aes(x = comm, fill = Phylum)) + geom_bar() + facet_grid(Gene~., scales = "free") +
  scale_fill_manual(values = as.character(ppal)) + xlab("Community") + ylab("Number of OTUs") +
  theme(panel.grid = element_blank())
  #scale_fill_manual(values = as.character(ppal), guide = guide_legend(reverse=TRUE))

ggsave(paste0("Community_composition_phylum_barplot.pdf"), 
       plot = p2, height = 21, width = 21, units = c("cm"), useDingbats = FALSE)

### Output summaries of community composition ###
dt <- data.table(d)
x <- dt[, .N, by=.(Gene, comm, Kingdom)]
write.table(x, file = "communities_by_kingdom.txt", sep = "\t", quote = F, col.names = NA)

x <- dt[, .N, by=.(Gene, comm, Kingdom, Phylum)]
write.table(x, file = "communities_by_phylum.txt", sep = "\t", quote = F, col.names = NA)
x1 <- dcast(x, Gene + Kingdom + Phylum ~ comm)

x <- dt[, .N, by=.(Gene, comm, Class)]
write.table(x, file = "communities_by_class.txt", sep = "\t", quote = F, col.names = NA)

}