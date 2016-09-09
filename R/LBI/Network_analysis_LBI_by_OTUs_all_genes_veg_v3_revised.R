### Network analysis by OTUs all genes plus veg ###

##### Spiec-easi example from
# https://github.com/zdk123/SpiecEasi

#library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(igraph)
library(qgraph)
library(phyloseq) # Use for filtering OTU data
library(data.table)

library(RColorBrewer)
colors <- brewer.pal(4, "Spectral")
pal <- colorRampPalette(colors) 
#devtools::install_github("hoesler/rwantshue")
#library(rwantshue)
#col_scheme <- iwanthue()

library("BiocParallel")
register(MulticoreParam(4))

library(gdata)
library(GGally)
library(Matrix)
library(network)
#library(rgexf)
#library(XML)

#source("H:/My Documents/PhD Research PFR folder/R stuff/ggnet2_mod.R")

###############################################################################

##### Analysis of LBI data #####

## Load taxonomic reference
# tax_ref <- read.table("D:/PhD_folder/Database stuff/New_taxonomy_from_PLOSONE_2015_fixed.txt", 
#                       header = TRUE, sep = "\t", quote = "", comment.char = "")
# tax_ref_2 <- read.table("D:/PhD_folder/Database stuff/New_taxonomy_from_PLOSONE_2015_hierarchy.txt", 
#                         header = TRUE, sep = "\t", quote = "", comment.char = "")
tax_ref <- read.table("H:/My Documents/PhD Research PFR folder/Database stuff/New_taxonomy_from_PLOSONE_2015_fixed.txt", 
                      header = TRUE, sep = "\t", quote = "", comment.char = "")
tax_ref_2 <- read.table("H:/My Documents/PhD Research PFR folder/Database stuff/New_taxonomy_from_PLOSONE_2015_hierarchy.txt", 
                        header = TRUE, sep = "\t", quote = "", comment.char = "")

# Remove quirks/questions in taxa ([= ...])
tax_ref$Order <- gsub("\\s[\\[\\]=()\\w\\s?.,]*", "", tax_ref$Order, perl = TRUE)
tax_ref$Class <- gsub("\\s[\\[\\]=()\\w\\s?.,]*", "", tax_ref$Class, perl = TRUE)
tax_ref$Phylum <- gsub("\\s[\\[\\]=()\\w\\s?.,]*", "", tax_ref$Phylum, perl = TRUE)
tax_ref$Kingdom <- gsub("\\s[\\[\\]=()\\w\\s?.,]*", "", tax_ref$Kingdom, perl = TRUE)
tax_ref_2$Order <- gsub("\\s[\\[\\]=()\\w\\s?.,]*", "", tax_ref_2$Order, perl = TRUE)
tax_ref_2$Class <- gsub("\\s[\\[\\]=()\\w\\s?.,]*", "", tax_ref_2$Class, perl = TRUE)
tax_ref_2$Phylum <- gsub("\\s[\\[\\]=()\\w\\s?.,]*", "", tax_ref_2$Phylum, perl = TRUE)
tax_ref_2$Superphylum <- gsub("\\s[\\[\\]=()\\w\\s?.,]*", "", tax_ref_2$Superphylum, perl = TRUE)
tax_ref_2$Subkingdom <- gsub("\\s[\\[\\]=()\\w\\s?.,]*", "", tax_ref_2$Subkingdom, perl = TRUE)
tax_ref_2$Kingdom <- gsub("\\s[\\[\\]=()\\w\\s?.,]*", "", tax_ref_2$Kingdom, perl = TRUE)

##### Load NCBI ranks (for subsetting data by taxon identifications)
NCBI_ranks_list <-  read.table("H:/My Documents/PhD Research PFR folder/Database stuff/NCBI_taxa_Orders_Genera_ranks.txt", 
                               header = TRUE, sep = "\t", quote = "", comment.char = "")
# NCBI_ranks_list <-  read.table("D:/PhD_folder/Database stuff/NCBI_taxa_Orders_Genera_ranks.txt", 
#                                header = TRUE, sep = "\t", quote = "", comment.char = "")

phyla <- unique(tax_ref_2$Phylum)
cclasses <- unique(tax_ref_2$Class)
subclasses <- unique(tax_ref_2$Subclass)
orders <- NCBI_ranks_list$Node[NCBI_ranks_list$Rank == "Order"]
#suborders <- NCBI_ranks_list$Node[NCBI_ranks_list$Rank == "Suborder"]
#superfams <- NCBI_ranks_list$Node[NCBI_ranks_list$Rank == "Superfamily"]
fams <- NCBI_ranks_list$Node[grepl("amily", NCBI_ranks_list$Rank)]
#subfams <- NCBI_ranks_list$Node[NCBI_ranks_list$Rank == "Subfamily"]
genera <- NCBI_ranks_list$Node[NCBI_ranks_list$Rank == "Genus"]

##### Load vegetation data #####
# setwd("D:/PhD_folder/LBI_miseq_analysis_stuff/LBI_veg_data/")
setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_veg_data/")
f <- Sys.glob("LBI_Trees_Saplings_SBA.txt")
veg.df <- read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = F)
veg.tax <- read.table("LBI_veg_taxonomy_list.txt", sep = "\t", header = TRUE, check.names = F) 

# Drop missing/empty samples
veg.df$CM30b <- NULL
veg.df$CM30b44 <- NULL
veg.df$CM30c58 <- NULL
veg.df <- veg.df[rowSums(veg.df) > 0, ] # Remove empty rows!!!!!

veg.df1 <- otu_table(veg.df, taxa_are_rows = T)
phy <- phyloseq(veg.df1)
# Remove OTUs not seen at least x times in at least y % of the samples???
# phy.sub = filter_taxa(phy, function(x) sum(x >= 1) > (0.2*length(x)), TRUE)
# phy.sub = filter_taxa(phy.sub, function(x) sum(x >= 1) < (0.8*length(x)), TRUE)
# ntaxa(phy.sub)
# veg.df <- as.data.frame(phy.sub)

# organise vegetation taxonomy
veg.df$Kingdom <- "Vegetation"
veg.df$Phylum <- veg.tax$Phylum[match(rownames(veg.df), veg.tax$Code)]
veg.df$Class <- veg.tax$Class[match(rownames(veg.df), veg.tax$Code)]
veg.df$Order <- veg.tax$Order[match(rownames(veg.df), veg.tax$Code)]
veg.df$Family <- veg.tax$Family[match(rownames(veg.df), veg.tax$Code)]
veg.df$Genus <- veg.tax$Species[match(rownames(veg.df), veg.tax$Code)]
veg.df$Gene <- paste("Vegetation")

##### Load miseq data #####
#setwd("D:/PhD_folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/")
setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/")
#setwd("C:/Documents and Settings/Andrew/Desktop/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/")
#setwd("J:/PhD/PhD Research/NZGL01401_analysis/LBI_U8_OTUs_analyses")

# envdata.bysubplot <- read.table("LBI_environmental_data/LBI_all_env_data_by_subplot.txt", 
#                                 sep = "\t", header = TRUE, row.names = 1)
# envdata.byplot <- read.table("LBI_environmental_data/LBI_all_env_data_by_plot.txt", 
#                              sep = "\t", header = TRUE, row.names = 1)

# normalisation function
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

# Analysis V3 - limited taxonomic groups per gene + veg
#rranks <- c("order", "family", "genus")
vals <- c("90")#,"85","80")
#genelist <- c("16S","18S","26S","ITS","ShCO1","FolCO1")
#genelist <- c("16S","18S","26S","ShCO1")
genelist <- c("16S","18S","ITS","ShCO1")
#genelist <- c("16S","18S")
outpath <- "Network_analysis/Networks_by_OTUs_all_genes_veg_v3_phy_20-80_16S-18SITSCOI/"

#for (rrank in rranks){
for (val in vals){
  
  outfile <- paste0(outpath, "LBI_all_genes_veg_ggnet_by_OTUs_summary_", val, ".txt")
  
  all.data <- data.frame()
  
  for (gene in genelist){
    
  f <- Sys.glob(paste0("OTUtables/", gene, "*otutable_min2_by_plot.txt"))
  df <- read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = F)

  #label<- sapply(strsplit(f, "norm/"), "[[", 2)
  label<- sapply(strsplit(f, "OTUtables/"), "[[", 2)
  genelabel <- sapply(strsplit(label, "_ME1.0"), "[[", 1)
  print(paste("starting", label, "..."))
  print(paste(gene, "OTUs before filtering: ", nrow(df)))
  cat(paste("\n\n", gene, "OTUs before filtering: ", nrow(df), "\n"), file = outfile, append = TRUE)
  #print(rrank)
  
  ### Filter dataframe in phyloseq to get common/abundant OTUs ###
  # (to reduce data to manageable size)
  df1 <- otu_table(df, taxa_are_rows = T)
  phy <- phyloseq(df1)
  # Remove OTUs not seen at least x times in at least y % of the samples
  if(gene == "ShCO1"|gene == "FolCO1"){
    #phy.sub = filter_taxa(phy, function(x) sum(x >= 4) >= (0.05*length(x)), TRUE)
    phy.sub = filter_taxa(phy, function(x) sum(x >= 1) > (0.2*length(x)), TRUE)
    phy.sub = filter_taxa(phy.sub, function(x) sum(x >= 1) < (0.8*length(x)), TRUE)
  }else{
    #phy.sub = filter_taxa(phy, function(x) sum(x >= 10) >= (0.05*length(x)), TRUE)
    phy.sub = filter_taxa(phy, function(x) sum(x >= 1) > (0.2*length(x)), TRUE)
    phy.sub = filter_taxa(phy.sub, function(x) sum(x >= 1) < (0.8*length(x)), TRUE)
  }
  ntaxa(phy.sub)
  df.sub <- as.data.frame(phy.sub)
  #df.sub <- df
  print(paste(gene, "OTUs after abundance filtering:", nrow(df.sub)))
  cat(paste(gene, "OTUs after abundance filtering:", nrow(df.sub), "\n"), file = outfile, append = TRUE)
 # df.sub <- df
   
  ### Get taxonomic data ###
  #match <- Sys.glob(paste0("Taxa_tables/", genelabel, "*s_nt_paths_new_taxonomy_table_2.txt"))
  #match <- Sys.glob(paste0("Taxa_tables/", g, "_taxonomy_table_new.txt"))
  if(gene == "16S"|gene == "18S"|gene == "26S"|gene == "ITS"){
    match <- Sys.glob(paste0("Taxa_tables/Various_BLAST_thresholds_percent/", gene, "_taxatable_", val, ".txt"))
    val <- sapply(strsplit(match, "table_|\\.txt"), "[[", 2)
  }else if(gene == "ShCO1"|gene == "FolCO1"){
    match <- Sys.glob(paste0("Taxa_tables/Various_BLAST_thresholds_percent/", gene, "_taxatable_70.txt"))
    #val <- sapply(strsplit(match, "table_|\\.txt"), "[[", 2)
  }
  print(paste("Threshold:", val))
  
  taxaTable <- read.table(match, sep="\t", header=T, row.names=1, quote = "", comment.char = "")
  taxa <- taxaTable[,c("Kingdom","Phylum","Class","Order","Family","Genus")]

  # Merge taxonomic data with OTU data
  df.tax <- merge(df.sub, taxa, by = "row.names")
  rownames(df.tax) <- df.tax$Row.names
  df.tax$Row.names <- NULL
  # Exclude OTUs without a useful identification 
  # (Only retain prokaryotes for 16S, eukaryotes for the other amplicons)  
  if(gene == "16S"){
    df.subset <- subset(df.tax, (grepl("BACTERIA|ARCHAEA", df.tax$Kingdom)))  
  }else if(gene == "18S"){
    #df.subset <- subset(df.tax, !(grepl("root|cellular organisms|No hits|Not assigned|BACTERIA|ARCHAEA", df.tax$Kingdom)))
    #df.subset <- subset(df.tax, (grepl("ANIMALIA|PROTOZOA|CHROMISTA|FUNGI", df.tax$Kingdom)))
    df.subset <- subset(df.tax, (grepl("PROTOZOA|CHROMISTA", df.tax$Kingdom)))
  }else if(gene == "26S"){
    df.subset <- subset(df.tax, (grepl("FUNGI", df.tax$Kingdom)))
    #df.subset <- subset(df.tax, (grepl("FUNGI|PROTOZOA|CHROMISTA", df.tax$Kingdom)))
  }else if(gene == "ITS"){
    df.subset <- subset(df.tax, (grepl("FUNGI", df.tax$Kingdom)))
  }else if(gene == "ShCO1"){
    df.subset <- subset(df.tax, (grepl("ANIMALIA", df.tax$Kingdom)))
    #df.subset <- subset(df.tax, (grepl("ANIMALIA|PROTOZOA|CHROMISTA|FUNGI", df.tax$Kingdom)))
  }else if(gene == "FolCO1"){
    df.subset <- subset(df.tax, (grepl("ANIMALIA", df.tax$Kingdom)))
    #df.subset <- subset(df.tax, (grepl("ANIMALIA|PROTOZOA|CHROMISTA", df.tax$Kingdom)))
  }
  # Drop OTUs with only Subkingdom or higher identification 
   df.subset <- df.subset[!(df.subset$Phylum %in% tax_ref_2$Superkingdom|
                              df.subset$Phylum %in% tax_ref_2$Kingdom|
                              df.subset$Phylum %in% tax_ref_2$Subkingdom), ]
  
#  df.subset <- df.subset[!(df.subset$Kingdom %in% tax_ref_2$Superkingdom), ]
 
  
  print(paste(gene, "OTUs after removing non-assigned and non-target taxa:", nrow(df.subset)))
  cat(paste(gene, "OTUs after removing non-assigned and non-target taxa:", nrow(df.subset), "\n"), file = outfile, append = TRUE)
      
  # Exclude probably bogus taxa
  df.subset <- subset(df.subset, !(grepl("Cnidaria|Brachiopoda|Echinodermata|Porifera", df.subset$Phylum)|
                                  grepl("Bivalvia|Teleostei|Elasmobranchii|Polyplacophora", df.subset$Class)|
                                  grepl("Nudibranchia|Crocodylia|Serpentes|Testudines|Carnivora|Gymnophiona|Lagomorpha|Rodentia|Serpentes|Scorpiones", 
                                          df.subset$Order)))
    
  print(paste(gene, "OTUs after removing bogus taxa:", nrow(df.subset)))
  cat(paste(gene, "OTUs after removing bogus taxa:", nrow(df.subset), "\n"), file = outfile, append = TRUE)
  
  ##### Drop non-Order/Family/Genus ids
#   if(rrank == "genus"){
#     ##### Network analysis by Genus
#     df.subset2 <- df.subset[df.subset$Genus %in% genera, ] 
#   }else if(rrank == "family"){
#     ##### Network analysis by Family
#     df.subset2 <- df.subset[df.subset$Family %in% fams, ] 
#   }else if(rrank == "order"){
#     ##### Network analysis by Order
#     df.subset2 <- df.subset[df.subset$Order %in% orders, ] 
#   }else if(rrank == "class"){
#     ##### Network analysis by Class
#     df.subset2 <- df.subset[df.subset$Class %in% cclasses, ] 
#   }       
#   
#  print(paste("OTUs after taxonomic rank filtering: ", nrow(df.subset2)))

  # Normalise datasets to a common scale (0-1)
  OTUtable <- df.subset[, 1:28]
  OTUtable.norm <- as.data.frame(sapply(OTUtable, normalize))
  df.subset.norm <- cbind(OTUtable.norm, df.subset[, 29:34])
  
  if(gene == "FolCO1"){
    gene <- "COI-650"
  }else if(gene == "ShCO1"){
    gene <- "COI-300"
  }
  
  df.subset.norm$Gene <- gene
    
  all.data <- rbind(all.data, df.subset.norm)
  
} # End of for gene in genelist


dim(all.data)
cat(paste("Total number of OTUs:", nrow(all.data), "\n"), file = outfile, append = TRUE)

# Replace any NAs
all.data[is.na(all.data)] <- 0

# Drop missing data columns
all.data$CM30b51 <- NULL
dim(veg.df)
dim(all.data)

# Combine miseq and veg data for taxonomy    
df.both <- rbind(all.data, veg.df)

# Normalise veg dataset to common scale (0-1)
vegtable <- veg.df[, grepl("CM|Plot|LB", colnames(veg.df))]
vegtable.norm <- as.data.frame(sapply(vegtable, normalize))
dim(vegtable.norm)

# Combine veg data with soil data for network inference
OTUtable.both <- rbind(all.data[, grepl("CM|Plot|LB", colnames(all.data))], vegtable.norm)
dim(OTUtable.both)
OTUtable.t <- t(OTUtable.both)


  ### Apply Spiec-Easi pipeline
  # Requires non-normalized count OTU table with samples on rows and OTUs in columns
  # Hence need to transpose input OTUtable
  se.mb <- spiec.easi(OTUtable.t, method='mb', lambda.min.ratio=5e-5, 
                      nlambda=25, icov.select.params=list(rep.num=50))
  #   se.gl <- spiec.easi(OTUtable.t, method='glasso', lambda.min.ratio=1e-2,
  #                       nlambda=25, icov.select.params=list(rep.num=50))
  #   sparcc <- sparcc(OTUtable.t)
    
  se.mb$opt.index
  se.mb$opt.sparsity
  #   se.gl$opt.index
  #   se.gl$opt.sparsity
  # 
  #   # Fix for occassional incorrect glasso result
  #   gl.note <- " "
  #   if(se.gl$opt.sparsity == 0){ # 1 = zero sparsity (no edges)
  #     se.gl$opt.index <- se.gl$opt.index + 1
  #     se.gl$opt.sparsity <- se.gl$sparsity[[se.gl$opt.index]]
  #     se.gl$refit <- wi2net(symmpart(se.gl$icov[[se.gl$opt.index]]))
  #     se.gl$refit[abs(se.gl$refit) > 0] <- 1
  #     gl.note <- "opt.sparsity = 0, opt.index increased by 1" 
  #     print(gl.note)
  #   }
    
    
    #############################################################################
    
    print("Creating Igraph objects...")
    ### Create unweighted igraph objects from path
    ig.mb <- graph.adjacency(se.mb$refit, mode='undirected')
    # ig.gl <- graph.adjacency(se.gl$refit, mode='undirected')
    
    #plot(sparcc$Cor)
  #   thr <- (max(sparcc$Cor[sparcc$Cor < 1]))/2 # Arbitrary threshold for sparcc graph
  #   sparcc.graph <- abs(sparcc$Cor) >= thr 
  #   ig.sparcc <- graph.adjacency(sparcc.graph, mode='undirected', diag=FALSE)
    
    ### Create weighted igraph objects from covariance magnitude data
    # mb$beta and gl$icov contain covariance/inverse covariance magnitude data
    # symmpart returns a symmetric matrix 
    # wi2net (qgraph) converts a precision matrix (inverse of covariance matrix/icov) to a partial correlation matrix
    # First, look at correlation values
    plot(symmpart(se.mb$beta[[se.mb$opt.index]]))
  #   plot(wi2net(symmpart(se.gl$icov[[se.gl$opt.index]])))
  #   plot(sparcc$Cor, xlim=c(-0.2, 0.5), ylim=c(-0.2, 0.5))
    # Then convert to igraph objects
    ig.mb.w <- graph.adjacency(symmpart(se.mb$beta[[se.mb$opt.index]]), mode = "undirected", weighted = TRUE, diag=FALSE)
  # ig.gl.w <- graph.adjacency(wi2net(symmpart(se.gl$icov[[se.gl$opt.index]])), mode = "undirected", weighted = TRUE, diag=FALSE)
  
    #plot(sparcc$Cor, xlim=c(-0.2, 0.5), ylim=c(-0.2, 0.5))
  #   sparcc.graph.w <- sparcc$Cor
  #   sparcc.graph.w[abs(sparcc.graph.w) < thr] <- 0 # Arbitrary threshold
  #   ig.sparcc.w <- graph.adjacency(sparcc.graph.w, weighted = TRUE, mode='undirected', diag=FALSE)
  

  
    # Set vertex attributes
    V(ig.mb.w)$OTU <- paste(rownames(df.both))
    V(ig.mb.w)$Gene <- paste(df.both$Gene)
    V(ig.mb.w)$Gene.K <- paste0(df.both$Gene, "_", df.both$Kingdom)
    V(ig.mb.w)$Genus <- paste(df.both$Genus)
    V(ig.mb.w)$Family <- paste(df.both$Family)
    V(ig.mb.w)$Order <- paste(df.both$Order)
    V(ig.mb.w)$Class <- paste(df.both$Class)
    V(ig.mb.w)$Phylum <- paste(df.both$Phylum)
    V(ig.mb.w)$Kingdom <- paste(df.both$Kingdom)
    V(ig.mb.w)$vsize <- rowMeans(clr(OTUtable.t, 1))+3 # sets size of vertex proportional to clr-mean?
    V(ig.mb.w)$degree <- igraph::degree(ig.mb.w)
  
    summary(ig.mb.w)
  
    # Plot using ggnet2!
    if(length(E(ig.mb.w)) == 0){
      print("No edges in graph!")
      cat(paste("No edges in graph!"), file = outfile, append = TRUE)
    }else{
      print("Plotting graphs...")
    
    #palette.phy = c("FUNGI" = "red", "FOO" = "green", "BAR" = "blue")
      
#     palette.gene = c("16S" = "#e46d1a", "18S" = "#377eb8",
#                        "26S" = "#4daf4a", "ITS" = "#984ea3", "COI-300" = "#f781bf", 
#                        "COI-650" = "#e41a1c", "Vegetation" = "#ffff33")
    palette.gene = c("BACTERIA" = "#e46d1a",
                     "ARCHAEA" = "#e46d1a",
                     "CHROMISTA" = "#377eb8",
                     "PROTOZOA" = "#377eb8",
                     "FUNGI" = "#4daf4a", 
                     "ANIMALIA" = "#f781bf", 
                     "Vegetation" = "#ffff33")  
          
    # MB graph
    net.mb <- as_data_frame(ig.mb.w, what = "both")
    #sign.mb <- sapply(E(ig.mb.w)$weight, function(x) ifelse(x < 0, 'black', 'grey'))
    lty.mb <- sapply(E(ig.mb.w)$weight, function(x) ifelse(x < 0, 2, 1))
    
    g.mb.gene <- ggnet2(ig.mb.w, node.color = V(ig.mb.w)$Kingdom, palette = palette.gene, node.alpha = 0.5, 
                     node.label = V(ig.mb.w)$Gene, node.size = V(ig.mb.w)$degree, 
                     label.size = 2, label.alpha = 0.5,
                     edge.size = abs(E(ig.mb.w)$weight*5), edge.alpha = 0.5, edge.lty = lty.mb,
                     #edge.color = sign.mb,
                     edge.color = c("color", "grey50"),
                     legend.position = "none")
    
    g.mb.g <- ggnet2(ig.mb.w, node.color = V(ig.mb.w)$Kingdom, palette = palette.gene, node.alpha = 0.5, 
                     node.label = V(ig.mb.w)$Genus, node.size = V(ig.mb.w)$degree, 
                     label.size = 2, label.alpha = 0.5,
                     edge.size = abs(E(ig.mb.w)$weight*5), edge.alpha = 0.5, edge.lty = lty.mb,
                     #edge.color = sign.mb,
                     edge.color = c("color", "grey50"),
                     legend.position = "none")
  
    g.mb.f <- ggnet2(ig.mb.w, node.color = V(ig.mb.w)$Kingdom, palette = palette.gene, node.alpha = 0.5, 
                     node.label = V(ig.mb.w)$Family, node.size = V(ig.mb.w)$degree, 
                     label.size = 2, label.alpha = 0.5,
                     edge.size = abs(E(ig.mb.w)$weight*5), edge.alpha = 0.5, edge.lty = lty.mb,
                     #edge.color = sign.mb,
                     edge.color = c("color", "grey50"),
                     legend.position = "none")
    
    g.mb.o <- ggnet2(ig.mb.w, node.color = V(ig.mb.w)$Kingdom, palette = palette.gene, node.alpha = 0.5, 
                     node.label = V(ig.mb.w)$Order, node.size = V(ig.mb.w)$degree, 
                     label.size = 2, label.alpha = 0.5,
                     edge.size = abs(E(ig.mb.w)$weight*5), edge.alpha = 0.5, edge.lty = lty.mb,
                     #edge.color = sign.mb,
                     edge.color = c("color", "grey50"),
                     legend.position = "none")
    
    g.mb.c <- ggnet2(ig.mb.w, node.color = V(ig.mb.w)$Kingdom, palette = palette.gene, node.alpha = 0.5, 
                     node.label = V(ig.mb.w)$Class, node.size = V(ig.mb.w)$degree, 
                     label.size = 2, label.alpha = 0.5,
                     edge.size = abs(E(ig.mb.w)$weight*5), edge.alpha = 0.5, edge.lty = lty.mb,
                     #edge.color = sign.mb,
                     edge.color = c("color", "grey50"),
                     legend.position = "none")
    
    g.mb.p <- ggnet2(ig.mb.w, node.color = V(ig.mb.w)$Kingdom, palette = palette.gene, node.alpha = 0.5, 
                     node.label = V(ig.mb.w)$Phylum, node.size = V(ig.mb.w)$degree, 
                     label.size = 2, label.alpha = 0.5,
                     edge.size = abs(E(ig.mb.w)$weight*5), edge.alpha = 0.5, edge.lty = lty.mb,
                     #edge.color = sign.mb,
                     edge.color = c("color", "grey50"),
                     legend.position = "none")
    
    g.mb.nl <- ggnet2(ig.mb.w, node.color = V(ig.mb.w)$Kingdom, palette = palette.gene, node.alpha = 0.5, 
                     node.size = V(ig.mb.w)$degree, 
                     label.size = 2, label.alpha = 0.5,
                     edge.size = abs(E(ig.mb.w)$weight*5), edge.alpha = 0.5, edge.lty = lty.mb,
                     #edge.color = sign.mb,
                     edge.color = c("color", "grey50"),
                     legend.position = "none")
    
  #   g.mb2 <- ggnet2(ig.mb.w, node.color = V(ig.mb.w)$Kingdom, palette = "Set1", node.alpha = 0.5, 
  #           node.label = V(ig.mb.w)$Genus, node.size = V(ig.mb.w)$degree, 
  #           label.size = 2.5, label.alpha = 0.5,
  #           edge.size = abs(E(ig.mb.w)$weight*5), edge.alpha = 0.5, edge.lty = lty.mb,
  #           #edge.color = sign.mb,
  #           edge.color = c("color", "grey50"),
  #           legend.position = "none")
  #   
  #   g.mb.nl <- ggnet2(ig.mb.w, node.color = V(ig.mb.w)$Kingdom, palette = "Set1", node.alpha = 0.5, 
  #                  node.size = round(V(ig.mb.w)$vsize, 1), 
  #                  label.size = 2.5, label.alpha = 0.5,
  #                  edge.size = abs(E(ig.mb.w)$weight*5), edge.alpha = 0.5, edge.lty = lty.mb,
  #                  #edge.color = sign.mb,
  #                  edge.color = c("color", "grey50"),
  #                  legend.position = "none")
  #   
  #   g.mb2.nl <- ggnet2(ig.mb.w, node.color = V(ig.mb.w)$Kingdom, palette = "Set1", node.alpha = 0.5, 
  #                   node.size = V(ig.mb.w)$degree, 
  #                   label.size = 2.5, label.alpha = 0.5,
  #                   edge.size = abs(E(ig.mb.w)$weight*5), edge.alpha = 0.5, edge.lty = lty.mb,
  #                   #edge.color = sign.mb,
  #                   edge.color = c("color", "grey50"),
  #                   legend.position = "none")
  
    # Print graphs
    pdf(file=paste0(outpath, "LBI_all_genes_veg_ggnet_by_OTUs_", val, "_v3_20-80.pdf"), width=25/2.5, height=25/2.5, useDingbats = FALSE)
    plot(g.mb.g)
    plot(g.mb.f)
    plot(g.mb.o)
    plot(g.mb.c)
    plot(g.mb.p)
    plot(g.mb.gene)
    plot(g.mb.nl)
    dev.off()
  
    ## Overall graph statistics
    print("Calculating graph statistics...")
    cat(paste("\n\n all genes network graph by OTUs", val, "all genes veg summary statistics\n"), file = outfile, append = TRUE)
    cat(paste("\nSparsity:"), file = outfile, append = TRUE)
    write.table(as.data.frame(se.mb$sparsity), sep = "\t", quote = FALSE, col.names = NA, file = outfile, append = TRUE)
    cat(paste("\nOpt.index:", se.mb$opt.index, "Opt.sparsity:", se.mb$opt.sparsity, ", Opt.lambda:", se.mb$opt.lambda), file = outfile, append = TRUE)
    cat(paste("\nVertices:", vcount(ig.mb.w), "Edges:", ecount(ig.mb.w), "\n"), file = outfile, append = TRUE)
    cat(paste("\nGraph density:", graph.density(ig.mb.w)), file = outfile, append = TRUE)
    cat(paste("\nTransitivity:", transitivity(ig.mb.w, type="global")), file = outfile, append = TRUE)
    cat(paste("\nNumber of clusters:", clusters(ig.mb.w)$no), file = outfile, append = TRUE)
    write.table(as.data.frame(clusters(ig.mb.w)$csize), sep = "\t", quote = FALSE, col.names = NA, file = outfile, append = TRUE)
    cat("\nDegree distribution:", file = outfile, append = TRUE)
    write.table(as.data.frame(degree.distribution(ig.mb.w)), sep = "\t", quote = FALSE, col.names = NA, file = outfile, append = TRUE)
  
    print("Calculating vertex and edge statistics...")  
    # Vertex statistics
    net.mb$vertices$degree <- igraph::degree(ig.mb.w)
    net.mb$vertices$betw <- igraph::betweenness(ig.mb)
    betw <- igraph::betweenness(ig.mb, normalized = TRUE)
    betw.norm <- (betw - min(betw))/(max(betw) - min(betw))
    net.mb$vertices$betw_norm <- betw.norm
  #   
    # Edge statistics
    net.mb$edges$betw <- igraph::edge.betweenness(ig.mb) # Using weighted graph gives error
    net.mb$edges$from.k <- net.mb$vertices$Kingdom[match(net.mb$edges$from, rownames(net.mb$vertices))]
    net.mb$edges$from.p <- net.mb$vertices$Phylum[match(net.mb$edges$from, rownames(net.mb$vertices))]
    net.mb$edges$from.c <- net.mb$vertices$Class[match(net.mb$edges$from, rownames(net.mb$vertices))]
    net.mb$edges$from.o <- net.mb$vertices$Order[match(net.mb$edges$from, rownames(net.mb$vertices))]
    net.mb$edges$from.f <- net.mb$vertices$Family[match(net.mb$edges$from, rownames(net.mb$vertices))]
    net.mb$edges$from.g <- net.mb$vertices$Genus[match(net.mb$edges$from, rownames(net.mb$vertices))] 
    net.mb$edges$from.gene <- net.mb$vertices$Gene[match(net.mb$edges$from, rownames(net.mb$vertices))] 
    net.mb$edges$to.k <- net.mb$vertices$Kingdom[match(net.mb$edges$to, rownames(net.mb$vertices))]
    net.mb$edges$to.p <- net.mb$vertices$Phylum[match(net.mb$edges$to, rownames(net.mb$vertices))]
    net.mb$edges$to.c <- net.mb$vertices$Class[match(net.mb$edges$to, rownames(net.mb$vertices))]
    net.mb$edges$to.o <- net.mb$vertices$Order[match(net.mb$edges$to, rownames(net.mb$vertices))]
    net.mb$edges$to.f <- net.mb$vertices$Family[match(net.mb$edges$to, rownames(net.mb$vertices))]
    net.mb$edges$to.g <- net.mb$vertices$Genus[match(net.mb$edges$to, rownames(net.mb$vertices))] 
    net.mb$edges$to.gene <- net.mb$vertices$Gene[match(net.mb$edges$to, rownames(net.mb$vertices))] 
    write.table(net.mb$vertices, file = paste0(outpath, "LBI_all_genes_veg_ggnet_by_OTUs_vertices_data_", val, "_v3_20-80.txt"), 
                sep = "\t", quote = FALSE, col.names = NA)
    write.table(net.mb$edges, file = paste0(outpath, "LBI_all_genes_veg_ggnet_by_OTUs_edges_data_", val, "_v3_20-80.txt"), 
                sep = "\t", quote = FALSE, col.names = NA)
  
    # Summary data plots
    pdf(file=paste0(outpath, "LBI_all_genes_veg_ggnet_by_OTUs_summary_plots_", val, "_v3_20-80.pdf"), width=15/2.5, height=10/2.5, useDingbats = FALSE)
  
    plot(degree.distribution(ig.mb.w), xlab = "node degree", ylab = "degree distribution", main = paste("all genes veg degree distribution, MB"))
    lines(degree.distribution(ig.mb.w))
    plot(igraph::betweenness(ig.mb, normalized = FALSE), ylab = "betweenness", xlab = "node", main = paste("all genes veg betweenness, MB")) 
    plot(igraph::edge.betweenness(ig.mb), ylab = "betweenness", xlab = "node", main = paste("all genes veg edge betweenness, MB")) 
  
    dev.off()
  
    }
    print(paste(val, "done!"))
    
} # End of for val in vals         
#} # End of for rrank in rranks

  