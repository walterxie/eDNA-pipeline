
library(RColorBrewer)
colors <- brewer.pal(4, "Spectral")
pal <- colorRampPalette(colors) 
#devtools::install_github("hoesler/rwantshue")
#library(rwantshue)
#col_scheme <- iwanthue()

library("BiocParallel")
register(MulticoreParam(4))

library(ggplot2)
library(reshape2)
library(phyloseq)

###############################################################################

env_data <- read.table("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_plot.txt", 
                       sep = "\t", header = TRUE, row.names = 1, check.names = F)
env_data_sort <- env_data[order(env_data$Elevation),]

##### Load vegetation data #####
#setwd("D:/PhD_folder/LBI_miseq_analysis_stuff/LBI_veg_data/")
setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_veg_data/")
f <- Sys.glob("LBI_Trees_Saplings_SBA.txt")
veg.df <- read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = F)
veg.tax <- read.table("LBI_veg_taxonomy_list.txt", sep = "\t", header = TRUE, check.names = F) 
rownames(veg.tax) <- veg.tax$Code

# Drop missing data columns
#df.subset2$CM30b51 <- NULL  
veg.df$CM30b <- NULL
veg.df$CM30b44 <- NULL
veg.df$CM30c58 <- NULL
veg.df <- veg.df[rowSums(veg.df) > 0, ]

veg.df <- otu_table(veg.df, taxa_are_rows = T)
samples <- sample_data(env_data_sort)
veg.tax <- tax_table(as.matrix(veg.tax))
phy.v <- phyloseq(veg.df, samples, veg.tax)

#veg.df$Species <- veg.tax$Species[match(rownames(veg.df), veg.tax$Code)]
#df.l <- melt(veg.df, variable.name = "Plot")
#df.l$Elevation <- env_data$Elevation[match(df.l$Plot, rownames(env_data))]

# Samples sorted by elevation, plants by name
plot_heatmap(na.omit(phy.v), 
             low="#66CCFF", high="#000033", na.value = "White", 
             taxa.label = "Species", sample.order = rownames(env_data_sort), 
             taxa.order = rev(rownames(veg.tax)))

# Samples sorted by elevation, plants by similarity 
plot_heatmap(na.omit(phy.v), 
             low="#66CCFF", high="#000033", na.value = "White", 
             taxa.label = "Species", sample.order = rownames(env_data_sort),
             axis.text.y = element_text(face="italic")) # italics don't work?

# Samples and plants both sorted by similarity 
plot_heatmap(na.omit(phy.v), method = "NMDS", distance = "bray", 
             low="#66CCFF", high="#000033", na.value = "White", 
             taxa.label = "Species", 
             taxa.order = rev(rownames(veg.tax)))
