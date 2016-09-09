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
library(RColorBrewer)
colors <- brewer.pal(3, "RdYlBu")
#colors <- brewer.pal(4, "Paired")
#colors <- brewer.pal(4, "Set3")
pal <- colorRampPalette(colors) 


#sourcepath <- "C:/Documents and Settings/Andrew/Desktop/R stuff/"
sourcepath <- "H:/My Documents/PhD Research PFR folder/R stuff/"
source(paste0(sourcepath, "beta1_function.R"))

#setwd("C:/Documents and Settings/Andrew/Desktop/LBI_miseq_analysis_stuff/")
setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/")

### Load environmental data ###
envdata.bysubplot <- read.table("LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_subplot.txt", 
                                sep = "\t", header = TRUE, row.names = 1)
envdata.byplot <- read.table("LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_plot.txt", 
                             sep = "\t", header = TRUE, row.names = 1)

### Combine subplot columns (for vegetation comparisons ###
combine.subplots <- function(df.list){
  print("Combining subplot columns...")
  n <- 1
  df.list2 <- list()
  for(df in df.list){
    colnames(df) <- sapply(strsplit(colnames(df), "-"), "[[", 1) # Strip subplot letter
    df2 <- data.frame(matrix(ncol = 0, nrow = nrow(df))) # Empty data.frame with required number of rows
    for(col in unique(colnames(df))){
      #print(col)
      cols <- df[grep(col, colnames(df))] # Find each pair of subplot columns
      cols1 <- as.data.frame(rowSums(cols))
      colnames(cols1) <- col
      df2 <- cbind(df2, cols1) 
    }
    df.list2[[n]] <- df2
    names(df.list2)[[n]] <- names(df.list)[[n]]
    n <- n + 1
  }
  return(df.list2)
}

### Subset tables based on counts of OTUs (to exclude bad PCR results) ###
check.counts <- function(df){
  thr <- 19 # Minimum OTU count threshold for sample inclusion
  df1 <- df
  df1[df1 > 0] <- 1 # Convert to presence/absence
  x <- colSums(df1) # OTU counts per sample
  df2 <- df[, x > thr] # Subset columns based on OTU counts
  return(df2)
}

### Check that row names match for mantel and procrustes comparisons ###
check.rows.match <- function(d1, d2){
  dlist <- list()
  d1.m <- as.matrix(d1)
  d2.m <- as.matrix(d2)
  d1.c <- d1.m[complete.cases(d1.m), ]
  d2.c <- d2.m[complete.cases(d2.m), ]
  if(length(rownames(d1.c)) != length(rownames(d2.c))){
    print("nrows different")
    keep <- intersect(rownames(d1.c), rownames(d2.c))
    d1.c <- as.dist(d1.c[keep, keep])
    d2.c <- as.dist(d2.c[keep, keep])
  } else {
    d1.c <- d1
    d2.c <- d2
  }
  dlist[[1]] <- d1.c
  dlist[[2]] <- d2.c
  return(dlist)
}

### Retrieves legend from plot ###
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

### Get OTU datasets to be compared ###
get.datasets <- function(genes.taxa){
  df.list <- list()
  n <- 1
  for(g.t in genes.taxa){
    g <- g.t[[1]]
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
    #for(taxon in taxa){
    #print(taxon)
    taxon <- g.t[[2]]
    if(taxon == "assigned"){
      dfs <- subset(df, !(grepl("root|cellular organisms|No hits|Not assigned", df$Kingdom)))
      }else if(taxon == "proks-euks"){
        if(g == "16S"){
          dfs <- subset(df, grepl("PROKARYOTA", df$Superkingdom))
        }else{
          dfs <- subset(df, grepl("EUKARYOTA", df$Superkingdom))
      }
    }else if(taxon == "bacteria"){
      dfs <- subset(df, grepl("BACTERIA", df$Kingdom))
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
    #label <- paste(g) # If taxon label not wanted
    label <- paste(g, taxon)
    print(label)
    names(df.list)[[n]] <- label
    n <- n + 1
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
      if((nrow(df) > 25) & (ncol(df) > 25)) {   # Ignore tables with few remaining samples/OTUs
        print(paste("dims:", dim(df)))
        df.t <- t(df)
        df.t <- df.t[order(rownames(df.t)), ] # Ensure samples in consistent order 
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

### Mantel comparisons ###
#do.mantel <- function(dist.list, metric){
do.mantel <- function(dist.list){
  n <- 1
  m.results <- list() # list for mantel results storage
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
    man <- mantel(d1, d2, permutations = 999)
    m.results[[n]] <- list(t1 = labels(x[[i]])[[1]], t2 = labels(x[[i]])[[2]], 
                           man_stat = man$statistic, man_sig = man$signif)
    n <- n + 1
  }
  return(m.results)
}

### Vegetation Mantel comparisons ###
do.veg.mantel <- function(dist.list, veg){
  n <- 1
  m.results <- list() # list for mantel results storage
  for(i in 1:length(dist.list)){
    print(labels(dist.list)[[i]])
    if(length(labels(veg) != length(labels(dist.list[[i]])))){ # Different numbers of samples
      print("Making samples match...")
      # Subset to shared samples
      keep <- intersect(labels(veg), labels(dist.list[[i]]))
      m1 <- as.matrix(veg)
      m2 <- as.matrix(dist.list[[i]])
      d1 <- as.dist(m1[keep, keep])
      d2 <- as.dist(m2[keep, keep])
    } else { # Same numbers of samples
      print("Samples match...")
      d1 <- veg
      d2 <- dist.list[[i]]
    }
    man <- mantel(d1, d2, permutations = 999)
    m.results[[n]] <- list(t1 = "vegetation", t2 = labels(dist.list)[[i]], 
                           man_stat = man$statistic, man_sig = man$signif)
    n <- n + 1
  }
  return(m.results)
}

### Procrustes comparisons ###
#do.procrustes <- function(dist.list, metric){
do.procrustes <- function(dist.list){  
  p.results <- list() # list to store procrustes results 
  pro.list <- list() # list to store procrustes objects (for plotting)
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
    proc <- procrustes(d1.mds, d2.mds, scale = TRUE, symmetric = TRUE)
    prot <- protest(d1.mds, d2.mds, scale = TRUE, symmetric = TRUE, permutations = 999)
    p.results[[n]] <- list(t1 = labels(x[[i]])[[1]], t2 = labels(x[[i]])[[2]], 
                         pro_ss = prot$ss, pro_corr = prot$t0, pro_sig = prot$signif) 
    pro.list[[n]] <- proc
    names(pro.list)[[n]] <- paste(labels(x[[i]])[[1]], "vs.", labels(x[[i]])[[2]])
    n <- n + 1
  }
  return(list(p.results, pro.list))
}

### Vegetation procrustes comparisons ###
do.veg.procrustes <- function(dist.list, veg){  
  p.results <- list() # list to store procrustes results 
  pro.list <- list() # list to store procrustes objects (for plotting)
  n <- 1 # counter for plot output
  for(i in 1:length(dist.list)){
    print(labels(dist.list)[[i]])
    if(length(labels(veg) != length(labels(dist.list[[i]])))){ # Different numbers of samples
      print("Making samples match...")
      # Subset to shared samples
      keep <- intersect(labels(veg), labels(dist.list[[i]]))
      m1 <- as.matrix(veg)
      m2 <- as.matrix(dist.list[[i]])
      d1 <- as.dist(m1[keep, keep])
      d2 <- as.dist(m2[keep, keep])
    } else { # Same numbers of samples
      print("Samples match...")
      d1 <- veg
      d2 <- dist.list[[i]]
    }
    # Generate MDS plots
    d1.mds <- metaMDS(d1)
    d2.mds <- metaMDS(d2)
    proc <- procrustes(d1.mds, d2.mds, scale = TRUE, symmetric = TRUE)
    prot <- protest(d1.mds, d2.mds, scale = TRUE, symmetric = TRUE, permutations = 999)
    p.results[[n]] <- list(t1 = "vegetation", t2 = labels(dist.list)[[i]], 
                           pro_ss = prot$ss, pro_corr = prot$t0, pro_sig = prot$signif) 
    pro.list[[n]] <- proc
    names(pro.list)[[n]] <- paste("vegetation", " vs. ", labels(dist.list)[[i]])
    n <- n + 1
  }
  return(list(p.results, pro.list))
}


### Generate procrustes plots ###
get.procrustes.plots <- function(pro.list){
  n <- 1
  plot.list <- list()
  for(i in 1:length(pro.list)){
    pro <- pro.list[[i]]
    pts <- data.frame(yMDS1 = pro$Yrot[,1], yMDS2 = pro$Yrot[,2], # rotated matrix i.e. d2.mds
                      xMDS1 = pro$X[,1], xMDS2 = pro$X[,2]) # target matrix i.e. d1.mds
    if(grepl("-", rownames(pts)[[1]])){
      pts <- merge(pts, envdata.bysubplot, by = "row.names")
    } else { # Subplots merged (vegetation comparisons)
      pts <- merge(pts, envdata.byplot, by = "row.names")
    }
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
      #ggtitle(paste0(letters[[i]], ". ", labels(x[[i]])[[1]], " vs. ", labels(x[[i]])[[2]])) #, ", ", metric, " distance")) #+
      ggtitle(paste0(letters[[i]], ". ", names(pro.list)[[i]]))
    #coord_fixed()
    legend <- get_legend(p)  # Get legend
    p <- p + theme(legend.position = "none")
    plot.list[[n]] <- ggplotGrob(p)
    n <- n + 1
  }
  return(list(plot.list, legend))
}

### Output procrustes plots as pdf ###
output.plots <- function(plots, file.label, metric){
  plot.list <- plots[[1]]
  legend <- plots[[2]]
  ### Split plot.list for printing if length > 6 (too many plots for one page) ###
  if (length(plot.list) <= 3) {  # up to 3 plots; 1 row
    pdf(file = paste0("LBI_", file.label, "_min2_norm_", metric, "_procrustes_compact_revised.pdf"),
        width = 20/2.54, height = 20*(1/3)/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=3, nrow=1))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  } else if ((length(plot.list) > 3) && (length(plot.list) <= 6)) {  # 4 to 6 plots; 2 rows
    pdf(file = paste0("LBI_", file.label, "_min2_norm_", metric, "_procrustes_compact_revised.pdf"),
        width = 20/2.54, height = 20*(2/3)/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=3, nrow=2))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  } else if ((length(plot.list) > 6) && (length(plot.list) <= 9)) {  # 7 to 9 plots; 3 rows
    pdf(file = paste0("LBI_", file.label, "_min2_norm_", metric, "_procrustes_compact_revised.pdf"),
        width = 20/2.54, height = 20/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=3, nrow=3))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  } else if ((length(plot.list) > 9) && (length(plot.list) <= 12)) {  # 10 to 12 plots; 4 rows
    pdf(file = paste0("LBI_", file.label, "_min2_norm_", metric, "_procrustes_compact_revised.pdf"),
        width = 20/2.54, height = 20*(1+1/3)/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=3, nrow=4))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  } else if ((length(plot.list) > 12) && (length(plot.list) <= 15)) {  # 12 to 15 plots; 5 rows
    pdf(file = paste0("LBI_", file.label, "_min2_norm_", metric, "_procrustes_compact_revised.pdf"),
        width = 20/2.54, height = 20*(1+2/3)/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=3, nrow=5))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  } else if (length(plot.list) > 15) {  # more than 15 plots; 4 cols, 6 rows (reduced size)
    pdf(file = paste0("LBI_", file.label, "_min2_norm_", metric, "_procrustes_compact_revised.pdf"),
        width = 20/2.54, height = 20*(1+2/3)/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=4, nrow=6))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  }
  dev.off()
} 

###############################################################################

metric <- "Jaccard"
#genes <- c("16S","18S","26S","ITS","ShCO1","FolCO1")
#taxa <- c("all OTUs","assigned","euks","proks","protists","fungi","animals","plants")
#taxa <- c("all OTUs")
#file.label <- "all_OTUs" # label for output file

### All genes comparisons
genes.taxa <- list(list("16S","all OTUs"),list("18S","all OTUs"),list("26S","all OTUs"),list("ITS","all OTUs"),list("ShCO1","all OTUs"),list("FolCO1","all OTUs"))
file.label <- "all_OTUs_all_genes"

### Different taxon within gene comparisons 
l <- list()
n <- 1
genes.taxa <- list(list("18S","fungi"),list("18S","protists"),list("18S","animals"),list("18S","plants"))
file.label <- "18S_x4"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1
genes.taxa <- list(list("26S","fungi"),list("26S","protists"),list("26S","animals"),list("26S","plants"))
file.label <- "26S_x4"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1
genes.taxa <- list(list("ShCO1","fungi"),list("ShCO1","protists"),list("ShCO1","animals"),list("ShCO1","plants"))
file.label <- "ShCO1_x4"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1
genes.taxa <- list(list("FolCO1","protists"),list("FolCO1","animals"),list("FolCO1","plants"),list("FolCO1","bacteria"))
file.label <- "FolCO1_x4"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1

### Same taxon between gene comparisons
genes.taxa <- list(list("18S","fungi"),list("26S","fungi"),list("ITS","fungi"),list("ShCO1","fungi"))
file.label <- "fungi_x4"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1
genes.taxa <- list(list("18S","fungi"),list("26S","fungi"),list("ITS","fungi"),list("ShCO1","fungi"),list("FolCO1","fungi"))
file.label <- "fungi_x5"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1
genes.taxa <- list(list("18S","protists"),list("26S","protists"),list("ShCO1","protists"),list("FolCO1","protists"))
file.label <- "protists_x4"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1
genes.taxa <- list(list("18S","protists"),list("26S","protists"),list("ITS","protists"),list("ShCO1","protists"),list("FolCO1","protists"))
file.label <- "protists_x5"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1
genes.taxa <- list(list("18S","animals"),list("26S","animals"),list("ShCO1","animals"),list("FolCO1","animals"))
file.label <- "animals_x4"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1
genes.taxa <- list(list("18S","animals"),list("26S","animals"),list("ITS","animals"),list("ShCO1","animals"),list("FolCO1","animals"))
file.label <- "animals_x5"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1
genes.taxa <- list(list("18S","plants"),list("26S","plants"),list("ShCO1","plants"),list("FolCO1","plants"))
file.label <- "plants_x4"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1
genes.taxa <- list(list("16S","bacteria"),list("FolCO1","bacteria"))
file.label <- "bacteria_x2"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1

### Selected between gene/taxon comparisons
genes.taxa <- list(list("16S","bacteria"),list("18S","protists"),list("18S","fungi"),list("18S","animals"),list("26S","fungi"),list("ITS","fungi"),list("ShCO1","animals"))
file.label <- "main_groups_genes_x7"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1
genes.taxa <- list(list("16S","bacteria"),list("18S","protists"),list("18S","fungi"),list("18S","animals"),list("26S","fungi"),list("ShCO1","animals"))
file.label <- "main_groups_genes_x6"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1
genes.taxa <- list(list("16S","bacteria"),list("18S","protists"),list("ITS","fungi"),list("ShCO1","animals"))
file.label <- "main_groups_genes_x4_v2"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1
genes.taxa <- list(list("16S","bacteria"),list("18S","protists"),list("18S","fungi"),list("18S","animals"))
file.label <- "main_groups_16S_18S"
l[[n]] <- list(genes.taxa, file.label)
n <- n + 1

genes.taxa <- list(list("16S","bacteria"),list("FolCO1","bacteria"),list("18S","fungi"),list("26S","fungi"),list("ITS","fungi"),list("ShCO1","fungi"),list("18S","protists"),list("26S", "protists"),list("ShCO1","protists"),list("FolCO1","protists"),list("18S","animals"),list("26S", "animals"), list("ShCO1","animals"),list("FolCO1","animals"),list("18S","plants"),list("26S", "plants"), list("ShCO1","plants"),list("FolCO1","plants"))
file.label <- "ALL_main_groups_all_genes"

for(z in l){
  genes.taxa <- z[[1]]
  file.label <- z[[2]]
  
df.list <- get.datasets(genes.taxa)
dist.list <- get.dissimilarity(df.list, metric)
m.result <- rbindlist(do.mantel(dist.list))
p.output <- do.procrustes(dist.list)
p.result <- rbindlist(p.output[[1]])
plots <- get.procrustes.plots(p.output[[2]])
output.plots(plots, metric, file.label)

write.table(m.result, 
            file = paste0("LBI_min2_norm_", file.label, "_mantel_results_revised.txt"), 
            sep = "\t", row.names = F, quote = F) 
write.table(p.result, 
            file = paste0("LBI_min2_norm_", file.label, "_procrustes_results_revised.txt"), 
            sep = "\t", row.names = F, quote = F) 

}
###############################################################################
### Vegetation comparisons ###
genes.taxa <- list(list("16S","all OTUs"),list("18S","all OTUs"),list("26S","all OTUs"),list("ITS","all OTUs"),list("ShCO1","all OTUs"),list("FolCO1","all OTUs"))
file.label <- "all_OTUs_all_genes_veg"

genes.taxa <- list(list("16S","bacteria"),list("18S","protists"),list("18S","fungi"),list("18S","animals"),list("26S","fungi"),list("ITS","fungi"),list("ShCO1","animals"),list("FolCO1","animals"))
file.label <- "main_groups_all_genes_veg"

genes.taxa <- list(list("16S","bacteria"),list("18S","protists"),list("18S","fungi"),list("18S","animals"),list("26S","fungi"),list("ShCO1","animals"))
file.label <- "main_groups_most_genes_veg"

genes.taxa <- list(list("16S","bacteria"),list("FolCO1","bacteria"),list("18S","fungi"),list("26S","fungi"),list("ITS","fungi"),list("ShCO1","fungi"),list("18S","protists"),list("26S", "protists"),list("ShCO1","protists"),list("FolCO1","protists"),list("18S","animals"),list("26S", "animals"), list("ShCO1","animals"),list("FolCO1","animals"),list("18S","plants"),list("26S", "plants"), list("ShCO1","plants"),list("FolCO1","plants"))
file.label <- "ALL_main_groups_all_genes_veg"

### Load vegetation data ###
f <- Sys.glob("LBI_veg_data/LBI_Trees_Saplings_SBA.txt")
veg <- read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
veg <- vegdist(t(veg), method = "jaccard", binary = TRUE)

df.list <- get.datasets(genes.taxa)
df.list <- combine.subplots(df.list) # Combine subplots for vegetation comparisons!
dist.list <- get.dissimilarity(df.list, metric)
m.result <- rbindlist(do.veg.mantel(dist.list, veg))
p.output <- do.veg.procrustes(dist.list, veg)
p.result <- rbindlist(p.output[[1]])
plots <- get.procrustes.plots(p.output[[2]])
output.plots(plots, metric, file.label)

write.table(m.result, 
            file = paste0("LBI_min2_", file.label, "_mantel_results_btw_genes_revised.txt"), 
            sep = "\t", row.names = F, quote = F) 
write.table(p.result, 
            file = paste0("LBI_min2_", file.label, "_procrustes_results_btw_genes_revised.txt"), 
            sep = "\t", row.names = F, quote = F) 

###############################################################################
### Heatmap of correlation values? 
d1 <- read.table("LBI_U8_OTUs_analyses/LBI_mantel_procrustes_revised/LBI_min2_norm_ALL_main_groups_all_genes_mantel_results_revised.txt", sep = "\t", header = T)
d1$measure <- "Mantel"
colnames(d1) <- gsub("man_", "", colnames(d1))
colnames(d1) <- gsub("stat", "corr", colnames(d1))
# d1$t1 <- factor(d1$t1, levels = d$t1, ordered = TRUE)
# d1$t2 <- factor(d1$t2, levels = rev(d$t2), ordered = TRUE)
#d <- d1

d2 <- read.table("LBI_U8_OTUs_analyses/LBI_mantel_procrustes_revised/LBI_min2_norm_ALL_main_groups_all_genes_procrustes_results_revised.txt", sep = "\t", header = T)
d2$measure <- "Procrustes"
d2$pro_ss <- NULL
colnames(d2) <- gsub("pro_", "", colnames(d2))
x <- d2$t1 # Flip t1 and t2
d2$t1 <- d2$t2
d2$t2 <- x
# d2$t1 <- factor(d2$t1, levels = d2$t1, ordered = TRUE)
# d2$t2 <- factor(d2$t2, levels = rev(d2$t2), ordered = TRUE)
#d <- d2

d <- rbind(d1,d2)

### Drop 18s and 26S plants?
#d <- d[!(grepl("18S plants|26S plants", d$t1)), ]
#d <- d[!(grepl("18S plants|26S plants", d$t2)), ]

### Taxonomic order
d$t1 <- factor(d$t1, levels = c("16S bacteria","COI-650 bacteria","18S fungi","26S fungi","ITS fungi",
                                "COI-300 fungi","18S protists","26S protists","COI-300 protists",
                                "COI-650 protists","18S animals","26S animals","COI-300 animals",
                                "COI-650 animals","18S plants","26S plants","COI-300 plants",
                                "COI-650 plants"), ordered = TRUE)
d$t2 <- factor(d$t2, levels = c("COI-650 plants","COI-300 plants","26S plants","18S plants","COI-650 animals",
                                "COI-300 animals","26S animals","18S animals","COI-650 protists","COI-300 protists",
                                "26S protists","18S protists","COI-300 fungi","ITS fungi","26S fungi","18S fungi",
                                "COI-650 bacteria","16S bacteria"), ordered = TRUE)
### Gene order
d$t1 <- factor(d$t1, levels = c("16S bacteria","18S fungi","18S protists","18S animals","18S plants","26S fungi",
                                "26S protists","26S animals","26S plants","ITS fungi","COI-300 fungi","COI-300 protists",
                                "COI-300 animals","COI-300 plants","COI-650 bacteria","COI-650 protists",
                                "COI-650 animals","COI-650 plants"), ordered = TRUE)
d$t2 <- factor(d$t2, levels = c("COI-650 plants","COI-650 animals","COI-650 protists","COI-650 bacteria","COI-300 plants",
                                "COI-300 animals","COI-300 protists","COI-300 fungi","ITS fungi","26S plants",
                                "26S animals","26S protists","26S fungi","18S plants","18S animals","18S protists",
                                "18S fungi","16S bacteria"), ordered = TRUE)
#d$t2 <- factor(d$t2, levels = rev(d$t2), ordered = TRUE)

### Cluster results? ###
z <- dcast(d1, t1 ~ t2, value.var = "corr")
rownames(z) <- z$t1
z$t1 <- NULL
d <- dist(as.matrix(z)) # Euclidean distance
h <- hclust(d)
plot(h)

d1$t1 <- factor(d1$t1, levels = h$labels[h$order], ordered = TRUE)
d1$t2 <- factor(d1$t2, levels = rev(h$labels[h$order]), ordered = TRUE)
ggplot(d1) + geom_tile(aes(x = t1, y = t2, fill = corr))



d$corr_sig <- paste0(round(d$corr, 2), "\n(", d$sig, ")")
d$corr_sig <- gsub("0.001", "<0.001", d$corr_sig)
d$sig_ind <- ifelse(d$sig == 0.001, "", d$sig)
d$sig_ind <- ifelse(d$sig_ind > 0.001 & d$sig_ind <= 0.01, "**", d$sig_ind)
d$sig_ind <- ifelse(d$sig_ind > 0.01 & d$sig_ind <= 0.05, "*", d$sig_ind)
d$sig_ind <- ifelse(d$sig_ind > 0.05, "ns", d$sig_ind)

ggplot(d) + geom_tile(aes(x = t1, y = t2, fill = corr)) +
  #geom_text(data = d, aes(x = t1, y = t2, label = d$corr_sig), size = 1.8) +
  geom_text(data = d, aes(x = t1, y = t2, label = d$sig_ind), size = 1.8) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7), 
        axis.text.y = element_text(size = 7), panel.grid = element_blank(), 
        legend.title = element_text(size = 7), legend.key.width = unit(0.5, "cm")) +
   scale_fill_gradient2(high = "#f46d43", mid = "#ffffbf", low = "#3288bd",
                         midpoint = 0.5, limit = c(0,1), name = "Correlation") +
#  scale_fill_gradient2(high = "#f1a340", mid = "#f7f7f7", low = "#998ec3",
#                       midpoint = 0.5, limit = c(0,1), name = "Correlation") +
#  scale_fill_gradient(high = "#f46d43", low = "#3288bd", limit = c(0,1), name = "Correlation") +
  xlab("") + ylab("")

### Get average values ###
d <- d[!(grepl("18S plants|26S plants",d$t1)|grepl("18S plants|26S plants",d$t2)), ]
d$g1 <- sapply(strsplit(as.character(d$t1), " "), "[[", 1)
d$g2 <- sapply(strsplit(as.character(d$t2), " "), "[[", 1)
d$gboth <- paste(d$g1, d$g2)
d$gboth <- plotSort(d$gboth)

d$tax1 <- sapply(strsplit(as.character(d$t1), " "), "[[", 2)
d$tax2 <- sapply(strsplit(as.character(d$t2), " "), "[[", 2)
d$taxboth <- paste(d$tax1, d$tax2)
d$taxboth <- plotSort(d$taxboth)

dt <- data.table(d)
dt[,mean(corr),by=.(measure, gboth)]
dt[,mean(corr),by=.(measure, taxboth)]

###########################################################################
d1 <- read.table("LBI_U8_OTUs_analyses/LBI_mantel_procrustes_revised/LBI_min2_norm_all_OTUs_all_genes_mantel_results_revised.txt", sep = "\t", header = T)
d1$measure <- "Mantel"
colnames(d1) <- gsub("man_", "", colnames(d1))
colnames(d1) <- gsub("stat", "corr", colnames(d1))

d2 <- read.table("LBI_U8_OTUs_analyses/LBI_mantel_procrustes_revised/LBI_min2_norm_all_OTUs_all_genes_procrustes_results_revised.txt", sep = "\t", header = T)
d2$measure <- "Procrustes"
d2$pro_ss <- NULL
colnames(d2) <- gsub("pro_", "", colnames(d2))
x <- d2$t1 # Flip t1 and t2
d2$t1 <- d2$t2
d2$t2 <- x

d <- rbind(d1,d2)

d$t1 <- gsub(" all OTUs", "", d$t1)
d$t2 <- gsub(" all OTUs", "", d$t2)

d$t1 <- factor(d$t1, levels = c("16S","18S","26S","ITS","COI-300","COI-650"), ordered = TRUE)
d$t2 <- factor(d$t2, levels = rev(c("16S","18S","26S","ITS","COI-300","COI-650")), ordered = TRUE)

ggplot(d) + geom_tile(aes(x = t1, y = t2, fill = corr)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7), 
        axis.text.y = element_text(size = 7), panel.grid = element_blank(), 
        legend.title = element_text(size = 7), legend.key.width = unit(0.5, "cm")) +
  scale_fill_gradient2(high = "#f46d43", mid = "#ffffbf", low = "#3288bd",
                       midpoint = 0.5, limit = c(0,1), name = "Correlation") +
  #  scale_fill_gradient2(high = "#f1a340", mid = "#f7f7f7", low = "#998ec3",
  #                       midpoint = 0.5, limit = c(0,1), name = "Correlation") +
  #  scale_fill_gradient(high = "#f46d43", low = "#3288bd", limit = c(0,1), name = "Correlation") +
  xlab("") + ylab("")

###########################################################################
### Vegetation heatmap ###
d1 <- read.table("LBI_U8_OTUs_analyses/LBI_veg_mantel_procrustes_revised/LBI_min2_ALL_main_groups_all_genes_veg_mantel_results_btw_genes_revised.txt", sep = "\t", header = T)
d1$measure <- "Mantel"
colnames(d1) <- gsub("man_", "", colnames(d1))
colnames(d1) <- gsub("stat", "corr", colnames(d1))

d2 <- read.table("LBI_U8_OTUs_analyses/LBI_veg_mantel_procrustes_revised/LBI_min2_ALL_main_groups_all_genes_veg_procrustes_results_btw_genes_revised.txt", sep = "\t", header = T)
d2$measure <- "Procrustes"
d2$pro_ss <- NULL
colnames(d2) <- gsub("pro_", "", colnames(d2))

d3 <- read.table("LBI_U8_OTUs_analyses/LBI_veg_mantel_procrustes_revised/LBI_min2_all_OTUs_all_genes_veg_mantel_results_btw_genes_revised.txt", sep = "\t", header = T)
d3$measure <- "Mantel"
colnames(d3) <- gsub("man_", "", colnames(d3))
colnames(d3) <- gsub("stat", "corr", colnames(d3))
d3$t2 <- paste(d3$t2, "all_OTUs")

d4 <- read.table("LBI_U8_OTUs_analyses/LBI_veg_mantel_procrustes_revised/LBI_min2_all_OTUs_all_genes_veg_procrustes_results_btw_genes_revised.txt", sep = "\t", header = T)
d4$measure <- "Procrustes"
d4$pro_ss <- NULL
colnames(d4) <- gsub("pro_", "", colnames(d4))
d4$t2 <- paste(d4$t2, "all_OTUs")

d <- rbind(d1, d2, d3, d4)

### Factors for vegetation heatmap
d$gene <- sapply(strsplit(as.character(d$t2), " "), "[[", 1)
d$taxon <- sapply(strsplit(as.character(d$t2), " "), "[[", 2)
d$gene <- factor(d$gene, levels = rev(c("16S","18S","26S","ITS","COI-300","COI-650")), ordered = TRUE)
d$taxon <- gsub("_", " ", d$taxon)
d$taxon <- factor(d$taxon, levels = c("all OTUs","bacteria","fungi","protists","animals","plants"), ordered = TRUE)

d$corr_sig <- paste0(round(d$corr, 2), "\n(", d$sig, ")")
d$corr_sig <- gsub("0.001", "<0.001", d$corr_sig)
d$sig_ind <- ifelse(d$sig == 0.001, "", d$sig)
d$sig_ind <- ifelse(d$sig_ind > 0.001 & d$sig_ind <= 0.01, "**", d$sig_ind)
d$sig_ind <- ifelse(d$sig_ind > 0.01 & d$sig_ind <= 0.05, "*", d$sig_ind)
d$sig_ind <- ifelse(d$sig_ind > 0.05, "ns", d$sig_ind)

ggplot(d, aes(x = taxon, y = gene, fill = corr)) + geom_tile() +
  facet_grid(. ~ measure) +
  #geom_text(data = d, aes(x = t1, y = t2, label = d$corr_sig), size = 1.8) +
  geom_text(data = d, aes(x = taxon, y = gene, label = d$sig_ind), size = 1.8) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7), 
        axis.text.y = element_text(size = 7), panel.grid = element_blank(), 
        legend.title = element_text(size = 7), legend.key.width = unit(0.5, "cm"),
        strip.background = element_rect(fill=NA, colour = NA)) +
  scale_fill_gradient2(high = "#f46d43", mid = "#ffffbf", low = "#3288bd",
                       midpoint = 0.5, limit = c(0,1), name = "Correlation") +
  #  scale_fill_gradient2(high = "#f1a340", mid = "#f7f7f7", low = "#998ec3",
  #                       midpoint = 0.5, limit = c(0,1), name = "Correlation") +
  #  scale_fill_gradient(high = "#f46d43", low = "#3288bd", limit = c(0,1), name = "Correlation") +
  xlab("") + ylab("") 

