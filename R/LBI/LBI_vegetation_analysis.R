
library(ggplot2)
library(reshape2)
library(scales)
library(data.table)
library(gridExtra)
library(plyr)
library(vegan)
library(vegetarian)

theme_set(theme_bw(base_size=9))

# Colour palette
library(RColorBrewer)
colors <- brewer.pal(4, "Spectral")
#colors <- brewer.pal(4, "Paired")
#colors <- brewer.pal(4, "Set3")
pal <- colorRampPalette(colors) 

### Retrieves legend from plot
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

path = "C:/Documents and Settings/Andrew/Desktop/"
path = "H:/My Documents/PhD Research PFR folder/"

source(paste0(path, "R stuff/beta1_function.R"))
source(paste0(path, "R stuff/multiplot_function.R"))

metadata <- read.table(paste0(path, "LBI_miseq_analysis_stuff/LBI_veg_data/LBI_veg_metadata.txt"), 
                       sep = "\t", header = TRUE, row.names = 1)
envdata.bysubplot <- read.table(paste0(path, "LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_subplot.txt"), 
                                sep = "\t", header = TRUE, row.names = 1)
envdata.byplot <- read.table(paste0(path, "LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_plot.txt"), 
                             sep = "\t", header = TRUE, row.names = 1)

#setwd("H:/My Documents/PhD Research PFR folder/")
setwd("C:/Documents and Settings/Andrew/Desktop/LBI_miseq_analysis_stuff/")
setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/")
#setwd("Y:/PhD/PhD Research/Sequencing_results/NZGL01401_analysis/LBI_miseq_analysis_stuff")

find_hull <- function(df) df[chull(df$MDS1, df$MDS2), ]

#fileNames <- Sys.glob("LBI_veg_data/*table.txt")
fileNames <- Sys.glob("LBI_veg_data/*SBA.txt")


### LBI vegetation MDS plots

for (f in fileNames) {
  
  df <- read.table(f, sep = "\t", header = TRUE, row.names = 1)  
  label<- sapply(strsplit(f, "/"), "[[", 2)
  #label <- sapply(strsplit(label, "_table"), "[[", 1)
  label <- sapply(strsplit(label, ".txt"), "[[", 1)
  print(paste("starting", label, "..."))
    
  # Exclude 3 plots not in Miseq data ?
  #df$CM30b <- NULL
  #df$CM30b44 <- NULL
  #df$CM30c58 <- NULL
  
  # Make distance matrices
  df.t <- t(df)
  dist.jac <- vegdist(df.t, method = "jaccard", binary = TRUE)
  dist.bray <- vegdist(df.t, method = "bray")
  dist.horn <- vegdist(df.t, method = "horn")
  dist.beta1 <- beta1(df.t, 1)
  attributes(dist.beta1)$Labels <- rownames(df.t)
  
  # Generate MDS plots
  mds.jac <- metaMDS(dist.jac)
  mds.bray <- metaMDS(dist.bray)
  mds.horn <- metaMDS(dist.horn)
  mds.beta1 <- metaMDS(dist.beta1)
  
  # Extract data for ggplot
  pts.jac <- as.data.frame(mds.jac$points)
  pts.bray <- as.data.frame(mds.bray$points)
  pts.horn <- as.data.frame(mds.horn$points)
  pts.beta1 <- as.data.frame(mds.beta1$points)
  
  str.jac <- mds.jac$stress
  str.bray <- mds.bray$stress
  str.horn <- mds.horn$stress
  str.beta1 <- mds.beta1$stress
  
  # Add metadata
  pts.jac <- merge(pts.jac, envdata.byplot, by = "row.names")
  pts.bray <- merge(pts.bray, envdata.byplot, by = "row.names")
  pts.horn <- merge(pts.horn, envdata.byplot, by = "row.names")
  pts.beta1 <- merge(pts.beta1, envdata.byplot, by = "row.names")
  
  pts.jac$Forest.code <- factor(pts.jac$Forest.code, levels = c("VS2", "VS3", "VS5", "WF7", "WF9", "WF11", 
                                                                "WF12", "WF13", "MF20", "x"), ordered = TRUE)
  pts.bray$Forest.code <- factor(pts.bray$Forest.code, levels = c("VS2", "VS3", "VS5", "WF7", "WF9", "WF11", 
                                                                "WF12", "WF13", "MF20", "x"), ordered = TRUE)
  pts.horn$Forest.code <- factor(pts.horn$Forest.code, levels = c("VS2", "VS3", "VS5", "WF7", "WF9", "WF11", 
                                                                "WF12", "WF13", "MF20", "x"), ordered = TRUE)
  pts.beta1$Forest.code <- factor(pts.beta1$Forest.code, levels = c("VS2", "VS3", "VS5", "WF7", "WF9", "WF11", 
                                                                "WF12", "WF13", "MF20", "x"), ordered = TRUE)
  
  # Get hulls for outlining sample groups?
  hulls.jac <- ddply(pts.jac, "Forest.Type", find_hull)
  hulls.bray <- ddply(pts.bray, "Forest.Type", find_hull)
  hulls.horn <- ddply(pts.horn, "Forest.Type", find_hull)
  hulls.beta1 <- ddply(pts.beta1, "Forest.Type", find_hull)
  
  plist <- list()
  
  p1 <- ggplot(pts.jac, aes(x = MDS1, y = MDS2, colour = Elevation)) + 
                  geom_point(aes(colour = Elevation, shape = Forest.code), size = 3, alpha = 1) +
                  #scale_shape_manual(values = c(5,4,2)) +
                  scale_shape_manual(values = c(15,16,17,0,1,2,5,6,3), name="Forest type") +
                  scale_colour_gradientn(colours = c("blue", "orange"), name="Elevation (m)") +
                  geom_text(aes(label = pts.jac$Row.names, colour = pts.jac$Elevation), size = 2, vjust = 2, alpha = 0.5) +
                  geom_polygon(data = hulls.jac, aes(mapping = hulls.jac$Forest.code, 
                                                     colour = hulls.jac$Elevation), alpha = 0.5, fill = NA) +
                  theme(panel.grid = element_blank(), plot.title = element_text(size = 8), legend.position = "none") +
                  xlab("") + ylab("") +
                  ggtitle(paste0("a. Forest composition, Jaccard distance (stress: ", paste(round(str.jac, 3)), ")"))
  
#   p2 <- ggplot(pts.bray, aes(x = MDS1, y = MDS2, colour = Elevation)) + 
#                   geom_point(aes(colour = Elevation, shape = Forest.code), size = 3, alpha = 1) +
#                   scale_shape_manual(values = c(15,16,17,0,1,2,5,6,3)) +
#                   geom_text(aes(label = pts.bray$Plot, colour = pts.bray$Elevation), size = 2, vjust = 2, alpha = 0.5) +
#                   geom_polygon(data = hulls.bray, aes(mapping = hulls.bray$Forest.Type, 
#                                                       colour = hulls.bray$Elevation), alpha = 0.5, fill = NA) +
#                   scale_colour_gradientn(colours = c("blue", "orange")) +
#                   theme(panel.grid = element_blank(), plot.title = element_text(size = 8), theme(legend.position = "none")) +
#                   xlab("") + ylab("") +
#                   ggtitle(paste0("b. ", label, ", Bray-Curtis (stress: ", paste(round(str.horn, 3)), ")"))
#                 
#   p3 <- ggplot(pts.horn, aes(x = MDS1, y = MDS2, colour = Elevation)) + 
#                   geom_point(aes(colour = Elevation, shape = Forest.code), size = 3, alpha = 1) +
#                   scale_shape_manual(values = c(15,16,17,0,1,2,5,6,3)) +
#                   geom_text(aes(label = pts.horn$Plot, colour = pts.horn$Elevation), size = 2, vjust = 2, alpha = 0.5) +
#                   geom_polygon(data = hulls.horn, aes(mapping = hulls.horn$Forest.Type, 
#                                                       colour = hulls.horn$Elevation), alpha = 0.5, fill = NA) +
#                   scale_colour_gradientn(colours = c("blue", "orange")) +
#                   theme(panel.grid = element_blank(), plot.title = element_text(size = 8), theme(legend.position = "none")) + 
#                   xlab("") + ylab("") +
#                   ggtitle(paste0("c. ", label, ", Morisita-Horn (stress: ", paste(round(str.horn, 3)), ")"))

  p2 <- ggplot(pts.beta1, aes(x = MDS1, y = MDS2, colour = Elevation)) + 
                  geom_point(aes(colour = Elevation, shape = Forest.code), size = 3, alpha = 1) +
                  scale_shape_manual(values = c(15,16,17,0,1,2,5,6,3), name="Forest type") +
                  geom_text(aes(label = pts.beta1$Row.names, colour = pts.beta1$Elevation), size = 2, vjust = 2, alpha = 0.5) +
                  geom_polygon(data = hulls.beta1, aes(mapping = hulls.beta1$Forest.code, 
                                                       colour = hulls.beta1$Elevation), alpha = 0.5, fill = NA) +
                  scale_colour_gradientn(colours = c("blue", "orange"),name="Elevation (m)") +
                  theme(panel.grid = element_blank(), plot.title = element_text(size = 8)) +
                  xlab("") + ylab("") +
                  ggtitle(paste0("b. Forest composition, 1-beta1 distance (stress: ", paste(round(str.beta1, 3)), ")"))
  
  legend <- get_legend(p2)  # Get legend
  
  p2 <- p2 + theme(legend.position = "none")
  
  plist <-list(p1, p2)

  pdf(file = paste0("LBI_veg_data/", label, "_MDS_plots_SBA_hulls_v3.pdf"),
      width = 22/2.54, height = 11/2.54, useDingbats = FALSE)
  args.list <- c(plist, list(ncol=2, nrow=1))
  print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  plist <- list() # reset plot list 
  dev.off()
  
  #multiplot(plist[[1]], plist[[2]], cols = 2)
  #dev.off()

}

