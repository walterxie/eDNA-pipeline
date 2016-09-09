library(ggplot2)
library(reshape2)
library(scales)
library(data.table)
library(gridExtra)
library(grid)
library(vegan)
library(vegetarian)

theme_set(theme_bw(base_size=8))

# Colour palette
library(RColorBrewer)
colors <- brewer.pal(4, "Spectral")
#colors <- brewer.pal(4, "Paired")
#colors <- brewer.pal(4, "Set3")
pal <- colorRampPalette(colors) 

#sourcepath = "J:/PhD/PhD Research/"
#sourcepath = "C:/DOcuments and Settings/Andrew/Desktop/"
sourcepath= "H:/My Documents/PhD Research PFR folder/"

source(paste0(sourcepath, "R stuff/beta1_function.R"))
#source(paste0(sourcepath, "R stuff/multiplot_function.R"))

find_hull <- function(df) df[chull(df$MDS1, df$MDS2), ]

setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/")
#setwd("C:/Documents and Settings/Andrew/Desktop/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/")
#setwd("J:/PhD/PhD Research/NZGL01401_analysis/LBI_U8_OTUs_analyses")

#fileNames <- Sys.glob("LBI_miseq_analysis_stuff/LBI_veg_data/*table.txt")

#tax_ref <- read.table("Silva_taxa_order_list.txt", header = TRUE, sep = "\t")

envdata.bysubplot <- read.table("LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_subplot.txt", 
                                sep = "\t", header = TRUE, row.names = 1)
envdata.byplot <- read.table("LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_plot.txt", 
                             sep = "\t", header = TRUE, row.names = 1)

#envdata.bysubplot$Forest.abbr <- sapply(strsplit(as.character(envdata.bysubplot$Forest.Type), ":"), "[[", 1)
#envdata.byplot$Forest.abbr <- sapply(strsplit(as.character(envdata.byplot$Forest.Type), ":"), "[[", 1)
#envdata.byplot$Plot <- rownames(envdata.byplot)

###############################################################################
### Functions for MDS plotting ###

### Retrieves legend from plot ###
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

### Get datasets for MDS plots ###
get.datasets <- function(genes.taxa){
  df.list <- list()
  n <- 1
  #for(g in genes){
  for(z in genes.taxa){
    g <- z[[1]]
    taxon <- z[[2]]
    print(paste(g, taxon))  
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
    #for(taxon in taxa){
      #print(taxon)
      if(taxon == "assigned"){
        dfs <- subset(df, !(grepl("root|cellular organisms|No hits|Not assigned", df$Kingdom)))
      #}else if(taxon == "proks-euks"){
      }else if(taxon == "prokaryotes"){
        dfs <- subset(df, grepl("PROKARYOTA", df$superkingdom))
      }else if(taxon == "eukaryotes"){
        dfs <- subset(df, grepl("EUKARYOTA", df$superkingdom))
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
      label <- paste(g, taxon)
      print(label)
      print(nrow(dfs))
      names(df.list)[[n]] <- label
      n <- n + 1
    #}
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

### Generate MDS plots ###
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
        scale_colour_gradientn(colours = c("blue", "orange"), limits = c(0, max(pts$Elevation))) +
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
output.plots <- function(plots.legend, file.label, metric){
  plot.list <- plots.legend[[1]]
  legend <- plots.legend[[2]]
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


###############################################################################
metric <- "Jaccard"

### Specify combinations of genes and taxa for MDS plotting ###
genes <- c("16S","18S","26S","ITS","ShCO1","FolCO1")
#taxa <- c("all OTUs","assigned","euks","proks","protists","fungi","animals","plants")
taxa <- c("all OTUs")
file.label <- "all_OTUs" # label for output file

genes <- c("18S")
taxa <- c("protists","fungi","animals","plants")
file.label <- "18S_subsets"

genes <- c("26S","ITS")
taxa <- c("fungi")
file.label <- "26S_ITS_fungi"

genes <- c("ShCO1","FolCO1")
taxa <- c("animals")
file.label <- "COI_animals"

genes <- c("18S","26S","ITS","ShCO1")
taxa <- c("fungi")
file.label <- "fungi_x4"

genes <- c("18S","26S","ShCO1","FolCO1")
taxa <- c("protists")
file.label <- "protists_x4"

genes <- c("18S","26S","ShCO1","FolCO1")
taxa <- c("animals")
file.label <- "animals_x4"

genes <- c("18S","26S","ShCO1","FolCO1")
taxa <- c("plants")
file.label <- "plants_x4"

#df.list <- get.datasets(genes, taxa)
genes.taxa <- list(list("16S","bacteria"),list("18S","protists"),list("18S","fungi"),
                   list("18S","animals"),list("26S","fungi"),list("ShCO1","animals"))
file.label <- "various_x6"
df.list <- get.datasets(genes.taxa)
dist.list <- get.dissimilarity(df.list, metric)
plots.legend <- get.mds.plots(dist.list, metric)
output.plots(plots.legend, file.label, metric)

