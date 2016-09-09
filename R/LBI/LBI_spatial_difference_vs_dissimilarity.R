library(data.table)
library(ggplot2)#
library(reshape2)#
library(scales)#
library(phyloseq)
#library(phytools) # For tree midpoint rooting
library(vegan)
library(vegetarian)
library(pez)
library(picante)

#library(doParallel)
#registerDoParallel(cores = 4)

theme_set(theme_bw(base_size = 8))

#setwd("C:/Documents and Settings/Andrew/Desktop/LBI_miseq_analysis_stuff")
#setwd("C:/Users/Andrew/PhD_folder/")
setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/")
#setwd("J:/PhD/PhD Research/NZGL01401_analysis/")

sourcepath = "H:/My Documents/PhD Research PFR folder/"
source(paste0(sourcepath, "R stuff/beta1_function.R"))

### Load environmnental data ###
envdata.bysubplot <- read.table("LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_subplot.txt", sep="\t", header=T, row.names=1)
envdata.byplot <- read.table("LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_plot.txt", 
                             sep = "\t", header = TRUE, row.names = 1)
envdata.bysubplot$Plot <- rownames(envdata.bysubplot)
envdata.byplot$Plot <- rownames(envdata.byplot)

### Path for phylogenetic trees ###
tree_path <- "LBI_U8_OTUs_analyses/LBI_U8_RAxML_trees_v2/FFT-NS2_min2_RAxML_trees/Prok_Euk_trees/"

### Calculate between-plot spatial distances ###
# http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Spherical Law of Cosines (slc)

# Converts degrees to radians
deg2rad <- function(deg) return(deg*pi/180)
# Get latitude/longitude coordinates from envdata
lat_long <- envdata.byplot[, c(5,6)]
lat_long.list <- as.list(as.data.frame(t(lat_long)))

# Distance calculations
R <- 6371000 # Earth mean radius [m]
x <- length(lat_long.list)
dlist <- list()
n <- 1
for(i in 1:(x-1)){
  for(j in 2:(x)){
    lat1 <- deg2rad(lat_long[i,][[1]])
    long1 <- deg2rad(lat_long[i,][[2]])
    lat2 <- deg2rad(lat_long[j,][[1]])
    long2 <- deg2rad(lat_long[j,][[2]])
    dist <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
    print(dist)
    res <- list("plot1" = rownames(lat_long)[[i]], "plot2" = rownames(lat_long)[[j]], 
                "lat1" = lat1, "lat2" = lat2, "long1" = long1, "long2" = long2, "dist" = dist)
    dlist[[n]] <- res
    n <- n + 1
  }
}

dists.btw <- rbindlist(dlist)
write.table(dists.btw, "LBI_U8_OTUs_analyses/within-between_plots/LBI_between-plot_distances.txt", 
            sep = "\t", quote = F, col.names = NA)

### sort pairs of plot names ###
plotSort <- function(x)
  sapply(lapply(strsplit(x, " "), sort), paste, collapse=" ")

### Load between plot distances ###
plot.dists <- read.table("LBI_U8_OTUs_analyses/Spatial_within_between_analysis/LBI_between-plot_distances.txt", 
                         sep = "\t", header = TRUE, row.names = 1)
plot.dists$pair <- paste(plot.dists$plot1, plot.dists$plot2)
plot.dists$pair <- plotSort(plot.dists$pair)

### Sort subplot pairs ###
strSort <- function(x)
  sapply(lapply(strsplit(x, NULL), sort), paste, collapse="")

### Load subplot distances ###
s <- read.table("LBI_U8_OTUs_analyses/Spatial_within_between_analysis/LBI_subplot_distances.txt", 
                sep = "\t", header = TRUE, row.names = 1)
subplot.dists <- data.frame(t(combn(names(s), 2)), dist=s[lower.tri(s)])
subplot.dists$pair <- paste0(subplot.dists$X1, subplot.dists$X2)
subplot.dists$pair <- strSort(subplot.dists$pair)

### Get OTU table data for each gene ###
get.data <- function(g){
  print(g)
  f <- Sys.glob(paste0("LBI_U8_OTUs_analyses/OTUtables/", g, "*otutable_min2.txt"))
  #f <- Sys.glob(paste0("LBI_U8_OTUs_analyses/OTUtables/", g, "*otutable_min2_by_plot.txt"))
  df <- read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  t <- Sys.glob(paste0("LBI_U8_OTUs_analyses/Taxa_tables/", g, "*s_nt_paths_new_taxonomy_table.txt"))
  taxonomy <- read.table(t, sep="\t", header=TRUE, row.names=1, quote = "", comment.char = "")
  df$Superkingdom <- taxonomy$Superkingdom[match(rownames(df), rownames(taxonomy))]
  if(g == "16S"){
    df <- subset(df, grepl("PROKARYOTA", df$Superkingdom))
  }else{
    df <- subset(df, grepl("EUKARYOTA", df$Superkingdom))
  }
  df$Superkingdom <- NULL
  print(dim(df))
  return(df)
}

### Calculate various multivariate distance metrics ###
get.dissimilarity <- function(df, metric){
  df.dist <- 0
  df.t <- t(df)
  print(paste("Calculating", metric))
  if(metric == "Jaccard"){
    df.dist <- vegdist(df.t, method = "jaccard", binary = TRUE)
    #df.dist <- vegdist(df.t, method = "jaccard", binary = FALSE)
  }else if(metric == "Bray-Curtis"){
    df.dist <- vegdist(df.t, method = "bray")
  }else if(metric == "Morisita-Horn"){
    df.dist <- vegdist(df.t, method = "horn")
  }else if(metric == "1-beta1"){
    df.dist <- beta1(df.t, 1)
    attributes(df.dist)$Labels <- rownames(df.t)
  }else if(grepl("unif", metric)){ # Unifrac distances
    if(g != "ShCO1"){ 
      t <- read.tree(Sys.glob(paste0(tree_path,"RAxML_bestTree.", g,"*_OTUs_MAFFT_rxGTR")))
    }else if(g == "ShCO1"){ # Different tree for ShCO1 eukaryotes
      t <- read.tree(Sys.glob(paste0(tree_path,"ShCO1_LBI_contigs_min2_eukaryotes_OTUs_MAFFT_FastTree.newick")))
    }
    phy <- phyloseq(otu_table(df.t, taxa_are_rows = FALSE), sample_data(envdata.bysubplot), phy_tree(t))
    #phy <- phyloseq(otu_table(df.t, taxa_are_rows = FALSE), sample_data(envdata.byplot), phy_tree(t))
    if(metric == "unwt-unif"){
      df.dist <- UniFrac(phy, weighted = FALSE, normalized = TRUE) 
    }else if(metric == "wt-unif"){
      df.dist <- UniFrac(phy, weighted = TRUE, normalized = TRUE)
    }
    #else if(g == "ShCO1"){
    #  m <- matrix(1, nrow = 2, ncol = 2) # Dummy data for ShCO1
    #  df.dist <- as.dist(m)
    #}
  }
  return(df.dist)
}

## Get all within-between pairwise distances ###
get.within.between.dists <- function(df.dist){
  dist.list <- melt(as.matrix(df.dist))
  dist.list$p1 <- gsub("-.", "", dist.list$Var1)
  dist.list$p2 <- gsub("-.", "", dist.list$Var2)
  dist.list$s1 <- sapply(strsplit(as.character(dist.list$Var1), "-"), "[[", 2)
  dist.list$s2 <- sapply(strsplit(as.character(dist.list$Var2), "-"), "[[", 2)
  dist.list$wb <- ifelse(dist.list$p1 == dist.list$p2, "within", "between")
  dist.list <- dist.list[dist.list$Var1 != dist.list$Var2, ] # Drop same sample distances
  return(dist.list)
}

### Get elevation differences between plots/subplots ###
get.elev.dist <- function(df.dist){
  dist.list <- melt(as.matrix(df.dist))
  #dist.list$e1 <- envdata.bysubplot$Elevation[match(dist.list$Var1, envdata.bysubplot$Plot)]
  #dist.list$e2 <- envdata.bysubplot$Elevation[match(dist.list$Var2, envdata.bysubplot$Plot)]
  dist.list$e1 <- envdata.byplot$Elevation[match(dist.list$Var1, envdata.byplot$Plot)]
  dist.list$e2 <- envdata.byplot$Elevation[match(dist.list$Var2, envdata.byplot$Plot)]
  dist.list$dist <- abs(dist.list$e1-dist.list$e2)
  dist.list <- dist.list[dist.list$Var1 != dist.list$Var2, ]
  return(dist.list)
}

### Get within-plot spatial distances ###
get.subplot.dists <- function(df.dist){
  dist.list <- melt(as.matrix(df.dist))
  dist.list$p1 <- gsub("-.", "", dist.list$Var1)
  dist.list$p2 <- gsub("-.", "", dist.list$Var2)
  dist.list$s1 <- sapply(strsplit(as.character(dist.list$Var1), "-"), "[[", 2)
  dist.list$s2 <- sapply(strsplit(as.character(dist.list$Var2), "-"), "[[", 2)
  dist.list <- dist.list[dist.list$p1 == dist.list$p2, ] # Drop between-plot distances
  dist.list <- dist.list[dist.list$s1 != dist.list$s2, ] # Drop same sample distances
  dist.list$spair <- strSort(paste0(dist.list$s1, dist.list$s2))
  dist.list$dist <- subplot.dists$dist[match(dist.list$spair, subplot.dists$pair)]
  return(dist.list)
}

### Get between-plot spatial distances ###
get.plot.dists <- function(df.dist){
  dist.list <- melt(as.matrix(df.dist))
  dist.list$p1 <- gsub("-.", "", dist.list$Var1)
  dist.list$p2 <- gsub("-.", "", dist.list$Var2)
  dist.list <- dist.list[dist.list$p1 != dist.list$p2, ] # Drop within-plot distances
  dist.list$pair <- plotSort(paste(dist.list$p1, dist.list$p2)) 
  dist.list$dist <- plot.dists$dist[match(dist.list$pair, plot.dists$pair)]
  return(dist.list)
}

### Get correlation statistics ###
get.corr <- function(dist.data){
  corlist <- list()
  i <- 1
  for(g in unique(dist.data$gene)){
    for(m in unique(dist.data$metric)){
      print(paste(g, m))
      x <- dist.data[(dist.data$gene==g & dist.data$metric==m),]
      if(nrow(x) > 10){
        t <- cor.test(x$value, x$dist, alternative="greater", method=c("pearson"))
        res <- list("metric" = m, "gene" = g, "t.stat" = t$statistic, "pval" = t$p.value,
                    "corr" = t$estimate[1], "df" = t$parameter, "95_lower" = t$conf.int[1], "95_upper" = t$conf.int[2], 
                    "alt" = t$alternative)
        corlist[[i]] <- res
        i <- i + 1
      }
    }
  }
  return(corlist)
}

###############################################################################

genes <- c("16S","18S","26S","ITS","ShCO1","FolCO1")
metrics <- c("Jaccard", "Bray-Curtis", "Morisita-Horn", "1-beta1", "wt-unif","unwt-unif")
#metrics <- list("wt-unif","unwt-unif")

dist.data <- data.frame()
for(g in genes){
  df <- get.data(g)
  for(metric in metrics){
    df.dist <- get.dissimilarity(df, metric)
    #d <- get.elev.dist(df.dist)
    d <- get.within.between.dists(df.dist)
    #d <- get.subplot.dists(df.dist)
    #d <- get.plot.dists(df.dist)
    if(nrow(d) > 0){
      d$gene <- g
      d$metric <- metric
    }
    dist.data <- rbind(dist.data, d)
  }
}
 
write.table(dist.data, file = "All_genes_elevation_vs_distance_by_plot.txt", sep = "\t", quote = FALSE, col.names = NA)
#write.table(dist.data, file = "All_genes_subplots_vs_distance.txt", sep = "\t", quote = FALSE, col.names = NA)
#write.table(dist.data, file = "All_genes_plots_vs_distance.txt", sep = "\t", quote = FALSE, col.names = NA)

corlist <- get.corr(dist.data)
cor.data <- rbindlist(corlist)
cor.data$cp1 <- paste0(round(cor.data$corr, 2), " (", round(cor.data$pval, 3), ")")
cor.data$cp <- gsub("\\(0\\)", "\\(<0.001\\)", cor.data$cp1)

write.table(cor.data, file = "All_genes_elevation_vs_distance_by_plot_corr_stats.txt", sep = "\t", quote = FALSE, col.names = NA)
#write.table(cor.data, file = "All_genes_subplots_vs_distance_corr_stats.txt", sep = "\t", quote = FALSE, col.names = NA)
#write.table(cor.data, file = "All_genes_plots_vs_distance_corr_stats.txt", sep = "\t", quote = FALSE, col.names = NA)

### Make some plots ###

#dist.data <- read.table("All_genes_elevation_vs_distance_by_subplot.txt", sep = "\t", header = TRUE, row.names = 1)

dist.data$gene <- gsub("ShCO1","COI-300", dist.data$gene)
dist.data$gene <- gsub("FolCO1","COI-650", dist.data$gene)
dist.data$gene <- factor(dist.data$gene, 
                         levels = c("16S","18S","26S","ITS","COI-300","COI-650"), ordered = TRUE)
dist.data$metric <- gsub("unwt-unif","unweighted Unifrac", dist.data$metric)
dist.data$metric <- gsub("wt-unif","weighted Unifrac", dist.data$metric)
dist.data$metric <- factor(dist.data$metric, 
                           levels = c("Jaccard","Bray-Curtis","Morisita-Horn","1-beta1","unweighted Unifrac","weighted Unifrac"),
                            ordered = TRUE)
dist.data$wb <- factor(dist.data$wb, levels = c("within", "between"), ordered = TRUE)

### Distance correlation plots ###
p <- ggplot(dist.data, aes(x = dist, y = value)) + 
  geom_point(shape = 1) + 
  geom_smooth(method = "lm", se = FALSE) + 
  ylab("Distance measure") + 
  #xlab("Elevation difference (m)") +
  xlab("Spatial distance (m)") +
  facet_grid(metric ~ gene, scales = "free") +
  theme(strip.background = element_blank(), panel.grid = element_blank())

### Add correlation data to plot? ###  
cor.data$gene <- gsub("ShCO1","COI-300", cor.data$gene)
cor.data$gene <- gsub("FolCO1","COI-650", cor.data$gene)
cor.data$gene <- factor(cor.data$gene, 
                         levels = c("16S","18S","26S","ITS","COI-300","COI-650"), ordered = TRUE)
cor.data$metric <- gsub("unwt-unif","unweighted Unifrac", cor.data$metric)
cor.data$metric <- gsub("wt-unif","weighted Unifrac", cor.data$metric)
cor.data$metric <- factor(cor.data$metric, 
                           levels = c("Jaccard","Bray-Curtis","Morisita-Horn","1-beta1","unweighted Unifrac","weighted Unifrac"),
                           ordered = TRUE)

#long cut way to find number of facets
len <- length(levels(dist.data$metric))*length(levels(dist.data$gene))
vars <- data.frame(expand.grid(levels(dist.data$metric), levels(dist.data$gene)))
colnames(vars) <- c("metric", "gene")
dat <- data.frame(x = rep(10, len), y = rep(0.5, len), vars) # x and y are coordinates relative to x/y axis scales
dat$labs <- cor.data$cp[match(interaction(cor.data$metric, cor.data$gene), interaction(dat$metric, dat$gene))]

p + geom_text(data = dat, aes(x, y, label=labs, group=NULL), size = 2)

### Within-between distances boxplot ###
ggplot(dist.data, aes(x = wb, y = value, color = wb)) + 
  geom_boxplot(outlier.shape = 1, outlier.colour = alpha("black", 0.25)) + 
  ylab("Distance measure") + xlab("") +
  facet_grid(metric ~ gene, scales = "free") +
  scale_colour_brewer(palette = "Set1") +
  theme(strip.background = element_blank(), panel.grid = element_blank(),
        legend.position = "none")




