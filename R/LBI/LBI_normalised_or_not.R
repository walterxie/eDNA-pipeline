library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(reshape2)
library(scales)
#library(phyloseq)
library(data.table)
library(vegan)
library(vegetarian)

theme_set(theme_bw(base_size=8))

#sourcepath <- "C:/Documents and Settings/Andrew/Desktop/R stuff/"
sourcepath <- "H:/My Documents/PhD Research PFR folder/R stuff/"
#sourcepath <- "J:/PhD/PhD Research/R stuff/"
source(paste0(sourcepath, "beta1_function.R"))

#setwd("H:/My Documents/PhD Research PFR folder/")
#setwd("C:/Documents and Settings/Andrew/Desktop/LBI_miseq_analysis_stuff/")
setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/")
#setwd("J:/PhD/PhD Research/NZGL01401_analysis/")

### Load environmental data ###
envdata.bysubplot <- read.table("LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_subplot.txt", 
                                sep = "\t", header = TRUE, row.names = 1)
envdata.byplot <- read.table("LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_plot.txt", 
                             sep = "\t", header = TRUE, row.names = 1)


### Retrieves legend from plot ###
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

### Get datasets (not normalised) ###
get.datasets <- function(genes, taxa){
  df.list <- list()
  n <- 1
  for(g in genes){
    #print(g)
    f <- Sys.glob(paste0("LBI_U8_OTUs_analyses/OTUtables/", g, "*otutable_min2.txt"))
    #f <- Sys.glob(paste0("LBI_U8_OTUs_analyses/OTUtables_norm/", g, "*otutable_min2_normalized.txt"))
    #f <- Sys.glob(paste0("LBI_U8_OTUs_analyses/OTUtables/", g, "*otutable_min2_by_plot.txt"))
    df <- read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    t <- Sys.glob(paste0("LBI_U8_OTUs_analyses/Taxa_tables/", g, "*s_nt_paths_new_taxonomy_table.txt"))
    taxonomy <- read.table(t, sep="\t", header=TRUE, row.names=1, quote = "", comment.char = "")
    df$Superkingdom <- taxonomy$Superkingdom[match(rownames(df), rownames(taxonomy))]
    df$Kingdom <- taxonomy$Kingdom[match(rownames(df), rownames(taxonomy))]
    if(g == "ShCO1"){ # Fix COI names
      g <- "COI-300"
    }
    if(g == "FolCO1"){
      g <- "COI-650"
    }
    for(taxon in taxa){
      #print(taxon)
      if(taxon == "assigned"){
        dfs <- subset(df, !(grepl("root|cellular organisms|No hits|Not assigned", df$Kingdom)))
      }else if(taxon == "proks-euks"){
        if(g == "16S"){
          dfs <- subset(df, grepl("PROKARYOTA", df$Superkingdom))
        }else{
          dfs <- subset(df, grepl("EUKARYOTA", df$Superkingdom))
        }
      }else if(taxon == "protists"){
        dfs <- subset(df, grepl("CHROMISTA|PROTOZOA", df$Kingdom))
      }else if(taxon == "fungi"){
        dfs <- subset(df, grepl("FUNGI", df$Kingdom))
      }else if(taxon == "animals"){
        dfs <- subset(df, grepl("ANIMALIA", df$Kingdom))
      }else if(taxon == "plants"){
        dfs <- subset(df, grepl("PLANTAE", df$Kingdom))
      }else if(taxon == "all OTUs"){
        dfs <- df # No subsetting
      }
      dfs$Superkingdom <- NULL
      dfs$Kingdom <- NULL
      df.list[[n]] <- dfs
      label <- paste(g, taxon)
      print(label)
      print(nrow(dfs))
      names(df.list)[[n]] <- label
      n <- n + 1
    }
  }
  return(df.list)
}

### Get datasets (normalised) ###
get.datasets.norm <- function(genes, taxa){
  df.list <- list()
  n <- 1
  for(g in genes){
    #print(g)
    #f <- Sys.glob(paste0("LBI_U8_OTUs_analyses/OTUtables/", g, "*otutable_min2.txt"))
    f <- Sys.glob(paste0("LBI_U8_OTUs_analyses/OTUtables_norm/", g, "*otutable_min2_normalized.txt"))
    #f <- Sys.glob(paste0("LBI_U8_OTUs_analyses/OTUtables/", g, "*otutable_min2_by_plot.txt"))
    df <- read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    t <- Sys.glob(paste0("LBI_U8_OTUs_analyses/Taxa_tables/", g, "*s_nt_paths_new_taxonomy_table.txt"))
    taxonomy <- read.table(t, sep="\t", header=TRUE, row.names=1, quote = "", comment.char = "")
    df$Superkingdom <- taxonomy$Superkingdom[match(rownames(df), rownames(taxonomy))]
    df$Kingdom <- taxonomy$Kingdom[match(rownames(df), rownames(taxonomy))]
    if(g == "ShCO1"){ # Fix COI names
      g <- "COI-300"
    }
    if(g == "FolCO1"){
      g <- "COI-650"
    }
    for(taxon in taxa){
      #print(taxon)
      if(taxon == "assigned"){
        dfs <- subset(df, !(grepl("root|cellular organisms|No hits|Not assigned", df$Kingdom)))
      }else if(taxon == "proks-euks"){
        if(g == "16S"){
          dfs <- subset(df, grepl("PROKARYOTA", df$Superkingdom))
        }else{
          dfs <- subset(df, grepl("EUKARYOTA", df$Superkingdom))
        }
      }else if(taxon == "protists"){
        dfs <- subset(df, grepl("CHROMISTA|PROTOZOA", df$Kingdom))
      }else if(taxon == "fungi"){
        dfs <- subset(df, grepl("FUNGI", df$Kingdom))
      }else if(taxon == "animals"){
        dfs <- subset(df, grepl("ANIMALIA", df$Kingdom))
      }else if(taxon == "plants"){
        dfs <- subset(df, grepl("PLANTAE", df$Kingdom))
      }else if(taxon == "all OTUs"){
        dfs <- df # No subsetting
      }
      dfs$Superkingdom <- NULL
      dfs$Kingdom <- NULL
      df.list[[n]] <- dfs
      label <- paste(g, taxon)
      print(label)
      print(nrow(dfs))
      names(df.list)[[n]] <- label
      n <- n + 1
    }
  }
  return(df.list)
}


### Calculate various multivariate distance metrics ###
get.dissimilarity <- function(df.list, metric){
  dist.list <- list()
  n <- 1
  for(df in df.list){
    if(nrow(df) > 0){  # Not an empty table
      print(summary(colSums(df)))
      df <- as.data.frame(df[, colSums(df != 0) > 5])  # Exclude any samples with fewer than x OTUs
      df <- as.data.frame(df[, colSums(df) > mean(colSums(df))*0.025])  # Exclude any samples with excessively low abundance
      df <- as.data.frame(df[rowSums(df) > 0, ]) # Exclude any empty rows    
      if((nrow(df) > 25) & (ncol(df) > 28)) {   # Ignore tables with few remaining samples/OTUs   
        print(paste("dims:", dim(df)))
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
        }
      }
      dist.list[[n]] <- df.dist
      names(dist.list)[[n]] <- names(df.list)[[n]]
      n <- n + 1
    }
  }
  return(dist.list)
}

### Get mds plots ###
get.mds.plots <- function(dist.list, metric){
  plot.list <- list()
  n <- 1
  for(dist in dist.list){  
    mds <- metaMDS(dist)
    pts <- as.data.frame(mds$points)
    str <- mds$stress
    ### Add metadata
    pts <- merge(pts, envdata.bysubplot, by = "row.names")
    #pts <- merge(pts, envdata.byplot, by = "row.names")
    #pts$Forest.code <- gsub("x", "Unknown", pts$Forest.code)
    pts$Forest.type <- factor(pts$Forest.code, levels = c("VS2", "VS3", "VS5", "WF7", "WF9", "WF11", 
                                                          "WF12", "WF13", "MF20", "Unknown"), ordered = TRUE)
    ### Get hulls for outlining sample groups?
    #hulls <- ddply(pts, "Forest.Type", find_hull)
    p <- ggplot(pts, aes(x = MDS1, y = MDS2)) + 
      geom_point(aes(colour = Elevation, shape = Forest.type), size = 3, alpha = 0.75) +
      #scale_shape_manual(values = c(5,4,2)) +
      scale_shape_manual(values = c(15,16,17,0,1,2,5,6,3,4)) +
      geom_text(aes(label = pts$Row.names, colour = pts$Elevation), size = 2, vjust = 2.5, alpha = 0.5) +
      geom_polygon(aes(mapping = pts$Plot, colour = pts$Elevation), alpha = 0.75) +
      #geom_polygon(data = hulls, aes(mapping = hulls$Forest.Type, colour = hulls$Elevation), fill = NA) +
      #scale_colour_gradientn(colours = pal(4)) +
      #scale_color_gradient(low = "#0000cc", high = "#ffff00") +
      scale_colour_gradientn(colours = c("blue", "orange")) +
      xlab("") + ylab("") +
      #scale_colour_gradientn(colours = c("#225ea8", "#41b6c4", "#a1dab4", "#ffff44")) +
      theme(panel.grid = element_blank(), plot.title = element_text(size = 8), 
            plot.margin = unit(c(0.1,0.25,0.25,0), "cm")) + labs(colour="Elevation (m)", shape="Forest type") +
      #scale_x_reverse() +
      ggtitle(paste0(letters[[n]], ". ", names(dist.list)[[n]], ", ", metric, " distance (stress: ", paste(round(str, 3)), ")"))
    
    legend <- get_legend(p)  # Get legend
    p <- p + theme(legend.position = "none")
    plot.list[[n]] <- ggplotGrob(p)
    names(plot.list)[[n]] <- names(dist.list)[[n]]
    n <- n + 1
  }
  return(list(plot.list, legend))
}

### Output MDS plots as pdf ###
output.mds.plots <- function(mds.plots, file.label, metric){
  plot.list <- mds.plots[[1]]
  legend <- mds.plots[[2]]
  ### Split plot.list for printing if length > 6 (too many plots for one page) ###
  if(length(plot.list) == 1){
    pdf(file = paste0("LBI_", file.label, "_min2_norm_", metric, "_MDS_plots_by_subplot_revised.pdf"),
        width = 11/2.54, height = 9/2.54, useDingbats = FALSE)
    print (grid.arrange(plot.list[[1]], legend, nrow = 1, ncol=2, widths=c(1, 0.2)))
  }else if(length(plot.list) == 2){
    pdf(file = paste0("LBI_", file.label, "_min2_norm_", metric, "_MDS_plots_by_subplot_revised.pdf"),
        width = 20/2.54, height = 9/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=2, nrow=1))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.11)))
  }else if(length(plot.list) > 2 && length(plot.list) < 5){
    pdf(file = paste0("LBI_", file.label, "_min2_norm_", metric, "_MDS_plots_by_subplot_revised.pdf"),
        width = 20/2.54, height = 18/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=2, nrow=2))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.11)))
  }else if(length(plot.list) > 4 && length(plot.list) < 7){
    pdf(file = paste0("LBI_", file.label, "_min2_norm_", metric, "_MDS_plots_by_subplot_revised.pdf"),
        width = 20/2.54, height = 27/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=2, nrow=3))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.11)))
  }else if(length(plot.list) > 6 && length(plot.list) < 9){
    pdf(file = paste0("LBI_", file.label, "_min2_norm_", metric, "_MDS_plots_by_subplot_revised.pdf"),
        width = 20/2.54, height = 36/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=2, nrow=4))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.11)))
  }
  dev.off()
} 

### Mantel comparisons ###
do.mantel <- function(dist.list){
  n <- 1
  m.results <- list() # list for mantel results storage
  #nn <- 1 # counter for plot output
  #nnn <- 1 # counter for results storage
  x <- combn(dist.list, 2, simplify = F) # Get all pairs of dist matrices
  for(i in 1:length(x)){
    print(labels(x[[i]]))
    if(length(labels(x[[i]][[1]])) != length(labels(x[[i]][[2]]))){ # Different numbers of samples
      print("Making samples match...")
      # Subset to shared samples
      keep <- intersect(labels(x[[i]][[1]]), labels(x[[i]][[2]]))
      m1 <- as.matrix(x[[i]][[1]])
      m2 <- as.matrix(x[[i]][[2]])
      d1 <- as.dist(m1[keep, keep])
      d2 <- as.dist(m2[keep, keep])
    } else { # Same numbers of samples
      print("Samples match...")
      d1 <- x[[i]][[1]]
      d2 <- x[[i]][[2]]
    }
    man <- mantel(d1, d2, permutations = 1999)
    m.results[[n]] <- list(t1 = labels(x[[i]])[[1]], t2 = labels(x[[i]])[[2]], 
                           man_stat = man$statistic, man_sig = man$signif)
    betadi
    n <- n + 1
  }
  return(m.results)
}


### Procrustes comparisons ###
do.procrustes <- function(dist.list){  
  p.results <- list() # list for procrustes results storage
  plot.list <- list() # list for plot storage
  n <- 1 # counter for plot output
  x <- combn(dist.list, 2, simplify = F) # Get all pairs of dist matrices
  for(i in 1:length(x)){
    print(labels(x[[i]]))
    if(length(labels(x[[i]][[1]])) != length(labels(x[[i]][[2]]))){ # Different numbers of samples
      print("Making samples match...")
      # Subset to shared samples
      keep <- intersect(labels(x[[i]][[1]]), labels(x[[i]][[2]]))
      m1 <- as.matrix(x[[i]][[1]])
      m2 <- as.matrix(x[[i]][[2]])
      d1 <- as.dist(m1[keep, keep])
      d2 <- as.dist(m2[keep, keep])
    } else { # Same numbers of samples
      print("Samples match...")
      d1 <- x[[i]][[1]]
      d2 <- x[[i]][[2]]
    }
    # Generate MDS plots
    d1.mds <- metaMDS(d1)
    d2.mds <- metaMDS(d2)
    pro <- protest(d1.mds, d2.mds, permutations = 1999)
    p.results[[n]] <- list(t1 = labels(x[[i]])[[1]], t2 = labels(x[[i]])[[2]], 
                           pro_ss = pro$ss, pro_corr = pro$t0, pro_sig = pro$signif) 
    pts <- data.frame(yMDS1 = pro$Yrot[,1], yMDS2 = pro$Yrot[,2], # rotated matrix i.e. d2.mds
                      xMDS1 = pro$X[,1], xMDS2 = pro$X[,2]) # target matrix i.e. d1.mds
    pts <- merge(pts, envdata.bysubplot, by = "row.names")
    r1 = acos(pro$rotation[1,1]) # X axis rotation (radians)
    r2 = r1 + (pi/2) # Y axis rotation (radians)
    p <- ggplot(pts) +
      geom_point(aes(x = yMDS1, y = yMDS2, colour = Elevation), shape = 1.5, size = 1.5, alpha = 0.75) + # rotated i.e. d2 (circles)
      geom_point(aes(x = xMDS1, y = xMDS2, colour = Elevation), shape = 2, size = 1, alpha = 0.75) + # target i.e. d1 (triangles)
      scale_shape(solid = FALSE) + xlab("") + ylab("") +
      geom_segment(aes(x = yMDS1, y = yMDS2, xend = xMDS1, yend = xMDS2, colour = Elevation), alpha = 0.75) +
      #arrow = arrow(length = unit(0.2,"cm")), alpha = 0.75) + 
      #geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") + 
      #geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
      #geom_abline(intercept = 0, slope = tan(r1), colour = "grey") + 
      #geom_abline(intercept = 0, slope = tan(r2), colour = "grey") +
      #geom_text(aes(x = xMDS1, y = xMDS2, label = pts.mds$Row.names, colour = Elevation), size = 1.5, vjust = 1.5, alpha = 0.5) + 
      scale_colour_gradientn(colours = c("blue", "orange")) + labs(colour="Elevation (m)") +
      theme(panel.grid = element_blank(), plot.title = element_text(size = 8), 
            plot.margin = unit(c(0.1,0.1,0.1,0), "cm"), legend.key.width = unit(0.65, "cm")) +
      ggtitle(paste0(letters[[i]], ". ", labels(x[[i]])[[1]], " vs. ", labels(x[[i]])[[2]])) #, ", ", metric, " distance")) #+
    #coord_fixed()
    
    legend <- get_legend(p)  # Get legend
    p <- p + theme(legend.position = "none")
    plot.list[[n]] <- ggplotGrob(p)
    names(plot.list)[[n]] <- paste0(labels(x[[i]])[[1]], "-", labels(x[[i]])[[2]])
    n <- n + 1
  }
  return(list(p.results, list(plot.list, legend)))
}

### Output procrustes plots as pdf ###
output.pro.plots <- function(pro.plots, file.label, metric){
  plot.list <- pro.plots[[1]]
  legend <- pro.plots[[2]]
  ### Split plot.list for printing if length > 6 (too many plots for one page) ###
  if (length(plot.list) < 7) {  # up to 6 plots; 2 rows
    pdf(file = paste0("LBI_", file.label, "_min2_", metric, "_procrustes.pdf"),
        width = 20/2.54, height = 20*(2/3)/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=3, nrow=2))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  } else if ((length(plot.list) > 6) && (length(plot.list) < 10)) {  # 7 to 9 plots; 3 rows
    pdf(file = paste0("LBI_", file.label, "_min2_", metric, "_procrustes.pdf"),
        width = 20/2.54, height = 20/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=3, nrow=3))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  } else if ((length(plot.list) > 9) && (length(plot.list) < 13)) {  # 10 to 12 plots; 4 rows
    pdf(file = paste0("LBI_", file.label, "_min2_", metric, "_procrustes.pdf"),
        width = 20/2.54, height = 20*(1+1/3)/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=3, nrow=4))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  } else if (length(plot.list) > 12) {  # more than 12 plots; 5 rows
    pdf(file = paste0("LBI_", file.label, "_min2_", metric, "_procrustes.pdf"),
        width = 20/2.54, height = 20*(1+2/3)/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=3, nrow=5))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  }
  dev.off()
} 

###############################################################################
### 
metric <- "Jaccard"
genes <- c("16S","18S","26S","ITS","ShCO1","FolCO1")
#taxa <- c("all OTUs","assigned","euks","proks","protists","fungi","animals","plants")
taxa <- c("all OTUs")
#file.label <- "18S" # label for output file

### Not normalised analysis ###
file.label <- "all_OTUs" # label for output file
df.list <- get.datasets(genes, taxa)

### Normalised analysis  ###
file.label <- "all_OTUs_norm" # label for output file
df.list.norm <- get.datasets.norm(genes, taxa)

dist.list <- get.dissimilarity(df.list, metric) # Distance matrices
mds.plots <- get.mds.plots(dist.list, metric) # MDS plots
output.mds.plots(mds.plots, metric, file.label) # Print MDS plots
m.result <- rbindlist(do.mantel(dist.list)) # Mantel tests
p.output <- do.procrustes(dist.list) # Procrustes tests
p.result <- rbindlist(p.output[[1]]) 
pro.plots <- p.output[[2]]
output.pro.plots(pro.plots, metric, file.label) # Print procrustes plots

write.table(m.result, file = paste0("LBI_U8_OTUs_analyses/LBI_", file.label, "_min2_", metric, "_mantel_results.txt"), 
            sep = "\t", row.names = F, quote = F) 
write.table(p.result, file = paste0("LBI_U8_OTUs_analyses/LBI_", file.label, "_min2_", metric, "_procrustes_results.txt"), 
            sep = "\t", row.names = F, quote = F) 

###############################################################################
df1 <- df.list[[1]]
df2 <- df.list[[2]]
df3 <- df.list.norm[[1]]
df4 <- df.list.norm[[2]]
#df1 <- df.list[[1]][1:20,]
#df2 <- df.list[[2]][1:10,]

dist1 <- vegdist(t(df1), method = "jaccard", binary = TRUE)
dist2 <- vegdist(t(df2), method = "jaccard", binary = TRUE)
dist3 <- vegdist(t(df3), method = "jaccard", binary = TRUE)
dist4 <- vegdist(t(df4), method = "jaccard", binary = TRUE)

# Mantel based on distance matrices
m12 <- mantel(dist1, dist2) # No match
m34 <- mantel(dist3, dist4) # Significant match
m13 <- mantel(dist1, dist3) # No match
m24 <- mantel(dist2, dist4) # No match

# Procrustes based on distance matrices
p12 <- protest(dist1, dist2) # No match
plot(p12)
p34 <- protest(dist3, dist4) # Significant match
plot(p34)
p13 <- protest(dist1, dist3) # Marginal match?
plot(p13)
p24 <- protest(dist2, dist4) # No match
plot(p24)

# MDS ordinations
mds1 <- metaMDS(dist1)
mds2 <- metaMDS(dist2)
mds3 <- metaMDS(dist3)
mds4 <- metaMDS(dist4)
# Normalised and non-normalised ordinations are identical
plot(mds1, type = "t")
plot(mds3, type = "t")
plot(mds2, type = "t")
plot(mds4, type = "t")

# Procrustes based on MDS ordinations
p.m12 <- protest(mds1, mds2) # No match
plot(p.m12)
p.m34 <- protest(mds3, mds4) # Significant match
plot(p.m34)
p.m13 <- protest(mds1, mds3) # Marginal match?
plot(p.m13)
p.m24 <- protest(mds2, mds4) # No match
plot(p.m24)

### Betadisper comparisons? ###
df1 <- df.list[[1]]
df2 <- df.list[[2]]
df3 <- df.list.norm[[3]]
df4 <- df.list.norm[[4]]
#df1 <- df.list[[1]][1:20,]
#df2 <- df.list[[2]][1:10,]

colnames(df1) <- paste0(colnames(df1), "_a")
colnames(df2) <- paste0(colnames(df2), "_b")
df1$id <- 1:nrow(df1)
df2$id <- 1:nrow(df2)
dfboth <- merge(df1, df2, all = TRUE)
rownames(df3) <- df3$id
dfboth$id <- NULL
dfboth[is.na(df3)] <- 0
g <- sapply(strsplit(colnames(dfboth), "_"), "[[", 2)
dist <- vegdist(t(dfboth), method = "jaccard", binary = TRUE)
b <- betadisper(dist, g)
plot(b)

