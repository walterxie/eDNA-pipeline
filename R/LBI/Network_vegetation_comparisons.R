
library(igraph)
library(GGally)
library(reshape2)
library(ggplot2)
library(data.table)
theme_set(theme_bw(base_size=8))

# Path to veg only network
setwd(paste0("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/Network_analysis/Networks_vegetation_only"))
e.file <- Sys.glob("LBI_veg_species_only_ggnet_MB_graphs_edges_data.txt")
v.file <- Sys.glob("LBI_veg_species_only_ggnet_MB_graphs_vertices_data.txt")

verts.mb.1 <- read.table(v.file, sep = "\t", header = TRUE)
edges.mb.1 <- read.table(e.file, sep = "\t", header = TRUE, row.names = 1)
#edges.mb$to.sp <- verts.mb$Species[match(edges.mb$to, verts.mb$X)]
#write.table(edges.mb, file = "LBI_veg_species_only_ggnet_MB_graphs_edges_data.txt", sep = "\t", quote = F, col.names = NA)

verts.mb.1$X <- verts.mb.1$Species
edges.mb.1$from <- edges.mb.1$from.sp
edges.mb.1$to <- edges.mb.1$to.sp
edges.mb.1$weight <- 1

g <- graph.data.frame(edges.mb.1, directed=FALSE, vertices=verts.mb.1)  

summary(g)
lty.mb <- sapply(E(g)$weight, function(x) ifelse(x < 0, 2, 1))
sign.mb <- sapply(E(g)$weight, function(x) ifelse(x < 0, "#377eb8","#ffff33"))
g.mb.gene <- ggnet2(g, node.color = "#ffff33", node.alpha = 0.5, 
                    node.label = V(g)$Species, node.size = V(g)$degree, 
                    label.size = 2, label.alpha = 0.5,
                    edge.size = abs(E(g)$weight*5), edge.alpha = 0.5, edge.lty = lty.mb,
                    edge.color = sign.mb,
                    #edge.color = c("color", "grey50"),
                    legend.position = "none")

plot(g.mb.gene)
glist <- list(g)

n <- 1
for(l in c("a","b","c","d")){
  # Path to other networks
  setwd(paste0("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/Network_analysis/Networks_by_OTUs_all_genes_veg_v3_phy_20-80_", l, "/"))
  e.file <- Sys.glob("*edges_data_90_v3_20-80.txt")
  v.file <- Sys.glob("*vertices_data_90_v3_20-80.txt")
  
  verts.mb <- read.table(v.file, sep = "\t", header = TRUE)
  edges.mb <- read.table(e.file, sep = "\t", header = TRUE, row.names = 1)
  verts.mb <- verts.mb[verts.mb$Gene == "Vegetation", ]
  edges.mb <- edges.mb[edges.mb$from.k == "Vegetation" & edges.mb$to.k == "Vegetation", ]
  verts.mb$X <- verts.mb$Genus
  edges.mb$from <- edges.mb$from.g
  edges.mb$to <- edges.mb$to.g
  #edges.mb$weight <- 1
  
  glist[[n]] <- graph.data.frame(edges.mb, directed=FALSE, vertices=verts.mb)
  n <- n + 1
}

summary(glist[[1]])
summary(glist[[2]])
summary(glist[[3]])
summary(glist[[4]])

g1 <- glist[[1]]
g2 <- glist[[3]]

sum(get.adjacency(g1) != get.adjacency(g2))/2
int <- graph.intersection(g1, g2)
ecount(g1)+ecount(g2)-2*ecount(int)

sum(get.adjacency(g) != get.adjacency(g2))/2
int <- graph.intersection(g, g2)
ecount(g)+ecount(g2)-2*ecount(int)


lty.mb <- sapply(E(g)$weight, function(x) ifelse(x < 0, 2, 1))
sign.mb <- sapply(E(g)$weight, function(x) ifelse(x < 0, "#377eb8","#ffff33"))
g.mb.gene <- ggnet2(g, node.color = V(g)$Gene, palette = palette.gene, node.alpha = 0.5, 
                    node.label = V(g)$Genus, node.size = V(g)$degree, 
                    label.size = 2, label.alpha = 0.5,
                    edge.size = abs(E(g)$weight*5), edge.alpha = 0.5, edge.lty = lty.mb,
                    edge.color = sign.mb,
                    #edge.color = c("color", "grey50"),
                    legend.position = "none")

plot(g.mb.gene)