library(ggplot2)
library(data.table)
library(igraph)
library(GGally)
#source("H:/My Documents/PhD Research PFR folder/R stuff/R_LBI_Miseq/IWantHue.R")
theme_set(theme_bw(base_size=8))

###############################################################################
# Compare composition/structure of networks/communities (requires OTU identifications to be retained between networks)
setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/Network_analysis/")

outfile <- "Network_analysis_20-80_community_consistency_v2.txt"

comms_a <- read.table("Networks_by_OTUs_all_genes_veg_v3_phy_20-80_a/LBI_all_genes_veg_ggnet_by_OTUs_fg_all_communities_vertices_data.txt", sep = "\t", header = T)
comms_b <- read.table("Networks_by_OTUs_all_genes_veg_v3_phy_20-80_c/LBI_all_genes_veg_ggnet_by_OTUs_fg_all_communities_vertices_data.txt", sep = "\t", header = T)
# comms_a <- read.table("Networks_by_OTUs_all_genes_veg_v3_phy_20-80_16S-18S26SCOI/LBI_all_genes_veg_ggnet_by_OTUs_fg_all_communities_vertices_data.txt", sep = "\t", header = T)
# comms_b <- read.table("Networks_by_OTUs_all_genes_veg_v3_phy_20-80_16S-18SITSCOI/LBI_all_genes_veg_ggnet_by_OTUs_fg_all_communities_vertices_data.txt", sep = "\t", header = T)

cat("\n\n\nNetwork community comparisons, 20-80 a. vs c.\n\n", file = outfile, append = TRUE)
#cat("\n\n\nNetwork community comparisons, 20-80 b. vs d.\n\n", file = outfile, append = TRUE)

### Drop fungal OTUs? ###
# comms_a <- comms_a[comms_a$Gene != "26S", ]
# comms_b <- comms_b[comms_b$Gene != "ITS", ]

### Need to change id column? ###
comms_a$X <- comms_a$OTU
comms_b$X <- comms_b$OTU

### Check how many communities in each dataset ###
max(comms_a$comm)
max(comms_b$comm)

### Calculate actual results ###
### Set up matrix to hold data ###
m <- matrix(nrow = max(comms_a$comm), ncol = max(comms_b$comm))
  
r_names <- list() # to hold rownames
c_names <- list() # to hold column names
  
for(i in 1:max(comms_a$comm)){
  comm_a <- comms_a[comms_a$comm == i, ]
  nodes_a <- unique(comm_a$X)
  r_names[i] <- paste0("comm.", i, ".", length(nodes_a))
  for (j in 1:max(comms_b$comm)) {
    comm_b <- comms_b[comms_b$comm == j, ]
    nodes_b <- unique(comm_b$X)
    c_names[j] <- paste0("comm.", j, ".", length(nodes_b))
    #rownames(m)$j <- paste(length(nodes_b)) 
    matches <- match(nodes_a, nodes_b)
    m[i,j] <- length(matches[!is.na(matches)])
    print(paste(i, "nodes_a:", length(nodes_a), j , "nodes_b:", length(nodes_b), "matching:", length(matches[!is.na(matches)])))
  }
}
colnames(m) <- c_names
rownames(m) <- r_names
write.table(m, file = outfile, sep = "\t", quote = F, col.names = NA, append = TRUE)
actual <- melt(m)

### Generate randomised results ###
results <- list()
for(n in 1:1000){
  print(paste("randomising", n))
  comms_a$comm <- sample(comms_a$comm)
  comms_b$comm <- sample(comms_b$comm)

  ### Set up matrix to hold data ###
  m <- matrix(nrow = max(comms_a$comm), ncol = max(comms_b$comm))
  
  r_names <- list() # to hold rownames
  c_names <- list() # to hold column names
  
  for(i in 1:max(comms_a$comm)){
    comm_a <- comms_a[comms_a$comm == i, ]
    nodes_a <- unique(comm_a$X)
    r_names[i] <- paste0("comm.", i, ".", length(nodes_a))
    for (j in 1:max(comms_b$comm)) {
      comm_b <- comms_b[comms_b$comm == j, ]
      nodes_b <- unique(comm_b$X)
      c_names[j] <- paste0("comm.", j, ".", length(nodes_b))
      #rownames(m)$j <- paste(length(nodes_b)) 
      matches <- match(nodes_a, nodes_b)
      m[i,j] <- length(matches[!is.na(matches)])
      #print(paste(i, "nodes_a:", length(nodes_a), j , "nodes_b:", length(nodes_b), "matching:", length(matches[!is.na(matches)])))
    }
  }
  colnames(m) <- c_names
  rownames(m) <- r_names
  
  results[[n]] <- melt(m)
}

res <- rbindlist(results)

res$type <- "predicted"
actual$type <- "actual"
all.data <- rbind(res, actual)
all.data$combo <- paste0(all.data$Var1, all.data$Var2)

get.color <- function(z){
  d <- all.data[grep(z, all.data$combo), ]
  print(paste(d$value[d$type == "actual"], max(d$value[d$type == "predicted"])))
  if(d$value[d$type == "actual"] > max(d$value[d$type == "predicted"])){
    #if(d$type[grep(max(d$value), d$value)] == "actual"){
    d$col <- "Actual"
  }else if(d$value[d$type == "actual"] == max(d$value[d$type == "predicted"])){
    d$col <- "Same"
  }else{
    d$col <- "Predicted"
  }
  #print(paste(d$col))
  return(d)
}  

new.data <- data.frame()
for(z in unique(all.data$combo)){
  d <- get.color(z)
  new.data <- rbind(new.data, d)
}

ggplot(new.data, aes(x = value, fill = col)) + geom_histogram(binwidth = 1) +
  #ggplot(res, aes(x = value)) + geom_freqpoly(binwidth = 1) +
  geom_vline(data = actual, aes(xintercept = value)) + 
  scale_fill_manual(values = c("#e41a1c","grey","#377eb8")) +
  facet_grid(Var1 ~ Var2, scales = "free") +
  theme(panel.grid = element_blank(), legend.position="none")

zzz <- all.data[all.data$Var2 == "comm.3.316",]
write.table(zzz, file = "blahblah.txt", sep = "\t", quote = F, col.names = NA)
write.table(new.data, file = "blah.txt", sep = "\t", quote = F, col.names = NA)

###############################################################################
### Compare networks with ITS vs 26S OTUs ###

verts_a <- read.table("Networks_by_OTUs_all_genes_veg_v3_phy_20-80_16S-18S26SCOI/LBI_all_genes_veg_ggnet_by_OTUs_vertices_data_90_v3_20-80.txt", sep = "\t", header = T)
edges_a <- read.table("Networks_by_OTUs_all_genes_veg_v3_phy_20-80_16S-18S26SCOI/LBI_all_genes_veg_ggnet_by_OTUs_edges_data_90_v3_20-80.txt", sep = "\t", header = T)
edges_a$X <- NULL
net_a <- graph.data.frame(edges_a, directed = FALSE, vertices = verts_a)

verts_b <- read.table("Networks_by_OTUs_all_genes_veg_v3_phy_20-80_16S-18SITSCOI/LBI_all_genes_veg_ggnet_by_OTUs_vertices_data_90_v3_20-80.txt", sep = "\t", header = T)
edges_b <- read.table("Networks_by_OTUs_all_genes_veg_v3_phy_20-80_16S-18SITSCOI/LBI_all_genes_veg_ggnet_by_OTUs_edges_data_90_v3_20-80.txt", sep = "\t", header = T)
edges_b$X <- NULL
net_b <- graph.data.frame(edges_b, directed = FALSE, vertices = verts_b)

### Drop fungal OTUs ###
verts_a <- verts_a[verts_a$Gene != "26S", ]
edges_a <- edges_a[edges_a$from.gene != "26S" & edges_a$to.gene != "26S", ]
net_a <- graph.data.frame(edges_a, directed=FALSE, vertices=verts_a)    

verts_b <- verts_b[verts_b$Gene != "ITS", ]
edges_b <- edges_b[edges_b$from.gene != "ITS" & edges_b$to.gene != "ITS", ]
net_b <- graph.data.frame(edges_b, directed=FALSE, vertices=verts_b)    

int <- intersection(net_a, net_b)
summary(int)


