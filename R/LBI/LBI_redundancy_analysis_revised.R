library(vegan)
library(ggplot2)
library(grid)
library(gridExtra)
library(VIF)
library(xtable)
library(scales)
library(corrgram)

source("H:/My Documents/PhD Research PFR folder/R stuff/vif_function.R")
#source("D:/PhD_folder/R stuff/vif_function.R")

theme_set(theme_bw(base_size=8))

setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/")
#setwd("D:/PhD_folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/")

###############################################################################
### Prepare environmental data for analysis ###

#envdata_bysubplot <- read.table("D:/PhD_folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_subplot.txt", sep="\t", header=T, row.names=1)
envdata_bysubplot <- read.table("LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_subplot.txt", sep="\t", header=T, row.names=1)

# Subset to variables and samples to be included in analysis
# drop CM30b51 and CM30b58, missing aspect data
# 10cm depth temp data shows better trend than surface, but missing a data point (LB1)
#envdata <- envdata_bysubplot[, c(4,5,8,9,14:22)]
envdata <- envdata_bysubplot[c(5:56), c(4,5,8,9,14:22)]
names(envdata)[names(envdata) == 'mean_T_surface'] <- 'Temp.'
names(envdata)[names(envdata) == 'Northness'] <- 'sin.aspect'
names(envdata)[names(envdata) == 'Eastness'] <- 'cos.aspect'

# Inspect data
#pdf("LBI_envdata_panel.pdf", width = 10, height = 8)
plot(envdata, gap = 0, panel = panel.smooth) 
#dev.off()
# Log transform chem variables
envdata[,c(5:8,9:11)] <- log(envdata[,c(5:8,9:11)], 2)
# Replace inf with zero
envdata[envdata == "-Inf"] <- 0
# Inspect data again; Note that certain variables are highly correlated (e.g. EC/Organic.C/Total.N) 
#pdf("LBI_envdata_panel_log2.pdf", width = 10, height = 8)
plot(envdata, gap = 0, panel = panel.smooth) 
#dev.off()

envdata.relabel <- envdata
envdata.relabel <- envdata.relabel[,c(1:4,13,5:12)]
colnames(envdata.relabel) <- gsub("\\.N", " N", colnames(envdata.relabel))
colnames(envdata.relabel) <- gsub("\\.P", " P", colnames(envdata.relabel))                                  
colnames(envdata.relabel) <- gsub("\\.C", " C", colnames(envdata.relabel))                                  
colnames(envdata.relabel) <- gsub("EC", "E.C.", colnames(envdata.relabel))
colnames(envdata.relabel) <- gsub("C N.", "C/N ", colnames(envdata.relabel))
plot(envdata.relabel, gap = 0, lower.panel = panel.smooth, upper.panel = panel.conf, cex.axis = 0.75, col.smooth = "purple", col = "#333333")
corrgram(envdata.relabel, gap = 0, lower.panel = panel.pts, upper.panel = panel.conf, cex.axis = 0.75, col = "#333333")

###############################################################################
### Functions for redundancy analysis ###

### Retrieves legend from plot
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

### Get datasets for MDS plots ###
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


### Redundancy analysis ###  
do.rda <- function(df.list, file.label){
  n <- 1
  #rdaLists <- list()
  rdaReducedList <- list()
  rdaForwardList <- list()
  rdaBackwardList <- list()
  #rdaReducedPlots <- list()
  #rdaForwardPlots <- list()
  #rdaBackwardPlots <- list()
  for(df in df.list){
    data.label <- names(df.list)[[n]]
    tableFile <- paste0("LBI_", file.label, "_RDA_tables.txt")
    logfile <- paste0("LBI_", file.label, "_RDA_log.txt")

    print(paste("OTU rows:", nrow(df)))
    cat(paste("\n\n", data.label, "\n"), file = logfile, append = T)
    cat(paste("\nNumber of rows and columns before filtering:", nrow(df), ncol(df), "\n"), file = logfile, append = T)
    print(dim(df))
    ### Drop samples not present in envdata 
    df$"CM30b51-B" <- NULL
    df$"CM30b51-M" <- NULL
    df$"CM30b58-C" <- NULL
    df$"CM30b58-I" <- NULL  
    
    if(nrow(df) > 0){  # Not an empty table
      #print(colSums(df))
      #cat(paste("\n", g, taxon, "colSums before filtering:\n"), file = logfile, append = T)
      #write.table(as.data.frame(colSums(df)), file = logfile, sep = "\t", quote = F, append = T)
      print(summary(colSums(df)))
      df <- as.data.frame(df[, colSums(df != 0) > 5])  # Exclude any samples with fewer than x OTUs
      df <- as.data.frame(df[, colSums(df) > mean(colSums(df))*0.025])  # Exclude any samples with excessively low abundance
      df <- as.data.frame(df[rowSums(df) > 0, ]) # Exclude any empty rows    
      if((nrow(df) > 25) & (ncol(df) > 28)) {   # Ignore tables with few remaining samples/OTUs   
      
        cat(paste("Number of rows and columns after filtering:", nrow(df), ncol(df), "\n"), file = logfile, append = T)
        print(dim(df))
        cat(paste("\n", file.label, "colSums after filtering:\n"), file = logfile, append = T)
        write.table(as.data.frame(colSums(df)), file = logfile, sep = "\t", quote = F, append = T)
        cat(paste("More than x rows: proceed with capscale", "\n\n"), file = logfile, append = T)
  
        communityMatrix <- t(df) # Transpose
        if((nrow(communityMatrix) != nrow(envdata))) { # Rows don't match
          both <- match(rownames(communityMatrix), rownames(envdata)) # Matching rows
          envdata1 <- envdata[both, ] # Subset to matching rows
        }else if((nrow(communityMatrix) == nrow(envdata))) { # Rows match
          envdata1 <- envdata # Continue
        }

        # Ensure data rows are in same order
        communityMatrix <- communityMatrix[order(rownames(communityMatrix)),]
        envdata1 <- envdata1[order(rownames(envdata1)),]
        
        ### Distance-based redundancy analysis using capscale (Code borrowed from Walter's pipeline)
        rda_table <- data.frame(row.names=c("Constrained","Unconstrained"))
        anova_table <- data.frame(row.names=colnames(envdata1))
        
        rda_0 <- capscale(communityMatrix ~ 1, envdata1, distance = "jaccard") # Empty model
        rda_1 <- capscale(communityMatrix ~ ., envdata1, distance = "jaccard") # Maximal model, use for auto model building
        print(head(summary(rda_1)))

        # Variance inflation factor - indicates highly correlated variables
        print(vif.cca(rda_1))
        
        rda_table$Inertia <- c(round(rda_1$CCA$tot.chi, 3), round(rda_1$CA$tot.chi, 3))
        rda_table$Proportion <- c(rda_1$CCA$tot.chi/rda_1$tot.chi, rda_1$CA$tot.chi/rda_1$tot.chi)
        rda_table$Proportion <- percent(rda_table$Proportion) # %
        
        constrained_inertia <- c()
        constrained_proportion <- c()
        # Test each variable individually
        for (i in 1:length(colnames(envdata1))) {
          rda_individual <- capscale(formula = as.formula(paste("communityMatrix", colnames(envdata1)[i], sep=" ~ ")), 
                                     envdata1, distance = "jaccard")
          constrained_inertia <- c(constrained_inertia, rda_individual$CCA$tot.chi)
          constrained_proportion <- c(constrained_proportion, rda_individual$CCA$tot.chi/rda_individual$tot.chi)
        }
        anova_table$Inertia <- constrained_inertia
        anova_table$Proportion <- percent(constrained_proportion)
        
        # Compute all the single terms in the scope argument that can be added or dropped from model, 
        # fit those models and compute a table of the changes in fit.
        add_1 <- add1(rda_0, scope=formula(rda_1), test="perm")
        anova_table$Pr <- add_1$Pr[-1] # 1st row is <none>
        colnames(anova_table)[3] <- "Pr($>$F)"
        
        # Build model after stepwise removal of collinear variables (vif >= 10; requires vif_function.R) 
        # variance inflation factor (VIF) quantifies the severity of multicollinearity in an ordinary least squares regression analysis. 
        env_reduced <- vif_func(in_frame = envdata1)
        print(env_reduced) # Remaining variables
        
        # Build model automatically from reduced variable set
        # (Unsure how to pass env_reduced variables to capscale formula; paste() doesn't work...)
        #	rda_reduced <- capscale(communityMatrix ~ slope.degree + Mean.Temp + Northness + Eastness + 
        #							pH + C.N.ratio + NO3.N + NH4.N + Olsen.P, envdata1, distance = "jaccard")
        rda_reduced <- capscale(formula = as.formula(paste("communityMatrix", paste(env_reduced, collapse=" + "), sep=" ~ ")), 
                                envdata1, distance = "jaccard")
        #if (verbose) head(summary(rda_reduced))
        head(summary(rda_reduced))
        anova_reduced <- anova(rda_reduced, by = "terms")
        
        rda_table$Inertia.R <- c(round(rda_reduced$CCA$tot.chi, 3), round(rda_reduced$CA$tot.chi, 3))
        rda_table$Proportion.R <- c(rda_reduced$CCA$tot.chi/rda_reduced$tot.chi, rda_reduced$CA$tot.chi/rda_reduced$tot.chi)
        rda_table$Proportion.R <- percent(rda_table$Proportion.R) # %
        
        anova_table$Reduced <- is.element(rownames(anova_table), rownames(anova_reduced))
        anova_table$Reduced[which(anova_table$Reduced==T)] <- anova_reduced$Pr[-length(anova_reduced$Pr)]
        colnames(anova_table)[4] <- "Reduced Pr($>$F)"
        
        # Choose a Model by Permutation Tests in Constrained Ordination using forward model selection
        rda_reduced_f <- ordistep(rda_0, scope = formula(rda_reduced), direction = "forward", permutations = 3999)
        rda_forward <- capscale(formula = as.formula(rda_reduced_f$call), data = envdata1, distance = "jaccard")
        #if (verbose) head(summary(rda_forward))
        head(summary(rda_forward))
        anova_forward <- anova(rda_forward, by = "terms")
        
        rda_table$Inertia.F <- c(round(rda_forward$CCA$tot.chi, 3), round(rda_forward$CA$tot.chi, 3))
        rda_table$Proportion.F <- c(rda_forward$CCA$tot.chi/rda_forward$tot.chi, rda_forward$CA$tot.chi/rda_forward$tot.chi)
        rda_table$Proportion.F <- percent(rda_table$Proportion.F) # %
        
        anova_table$Forward <- is.element(rownames(anova_table), rownames(anova_forward))
        anova_table$Forward[which(anova_table$Forward==T)] <- anova_forward$Pr[-length(anova_forward$Pr)]
        colnames(anova_table)[5] <- "Forward Pr($>$F)"
        
        # Choose a Model by Permutation Tests in Constrained Ordination using backward model selection
        rda_reduced_b <- ordistep(rda_reduced, scope = formula(rda_0), direction = "backward", permutations = 3999)
        rda_backward <- capscale(formula = as.formula(rda_reduced_b$call), data = envdata1, distance = "jaccard")
        #if (verbose) head(summary(rda_backward))
        head(summary(rda_backward))
        anova_backward <- anova(rda_backward, by = "terms")
        
        rda_table$Inertia.B <- c(round(rda_backward$CCA$tot.chi, 3), round(rda_backward$CA$tot.chi, 3))
        rda_table$Proportion.B <- c(rda_backward$CCA$tot.chi/rda_backward$tot.chi, rda_backward$CA$tot.chi/rda_backward$tot.chi)
        rda_table$Proportion.B <- percent(rda_table$Proportion.B) # %
        
        anova_table$Backward <- is.element(rownames(anova_table), rownames(anova_backward))
        anova_table$Backward[which(anova_table$Backward==T)] <- anova_backward$Pr[-length(anova_backward$Pr)]
        colnames(anova_table)[6] <- "Backward Pr($>$F)"
        
        anova_table[anova_table == 0] <- ""
        anova_table$Proportion <- gsub("%", "\\\\%", anova_table$Proportion)
        
        print(xtable(anova_table, caption = paste("Distance-based redundancy analysis ANOVA tests for", data.label), 
                     label = paste("tab:rdaAnova", data.label, sep = "_"), caption.placement = "top"), 
              sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)
        
        rda_table$Proportion <- gsub("%", "\\\\%", rda_table$Proportion)
        rda_table$Proportion.R <- gsub("%", "\\\\%", rda_table$Proportion.R)
        rda_table$Proportion.F <- gsub("%", "\\\\%", rda_table$Proportion.F)
        rda_table$Proportion.B <- gsub("%", "\\\\%", rda_table$Proportion.B)
        
        print(xtable(rda_table, caption = paste("Distance-based redundancy analysis of environmental variables and", data.label), 
                     label = paste("tab:rda", data.label, sep = "_"), caption.placement = "top"), 
              sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)
        
        rdaReducedList[[ data.label ]] <- rda_reduced
        rdaForwardList[[ data.label ]] <- rda_forward
        rdaBackwardList[[ data.label ]] <- rda_backward
        #rdaList <- list("Reduced" = rdaReduced, "Forward" = rdaForward, "Backward" = rdaBackward)
      }
    }
    n <- n + 1
  }  
  rda.lists <- list("Reduced" = rdaReducedList, "Forward" = rdaForwardList, "Backward" = rdaBackwardList)
  #rdaLists[[n]] <- rdaList
  #names(rdaLists)[[n]] <- names(df.list)[[n]]
  return(rda.lists)
}    

### Make biplots using ggplot ###
# sp = species scores, wa = site scores, bp = biplot arrows, lc = linear constraints 
get.rda.plots <- function(rda.pick, file.label){ 
  plot.list <- list()
  n <- 1
  for(rda in rda.pick){
    #lab <- names(rdaList)[[ n ]]
    #j <- 1
    #for(rda in rdaList){
      #taxon <- labels(rdaList)[[ j ]]
      sites <- as.data.frame(scores(rda, display = "wa"))
      sites$Elevation <- envdata$Elevation[match(rownames(sites), rownames(envdata))]
      sites$Plot <- gsub("-[A-Z]", "", rownames(sites)) 
      sites$shortIDs <- gsub("(CM30|CM31|Plot)", "", rownames(sites))
      biplots <- as.data.frame(scores(rda, display = "bp", scaling = "symmetric"))
      biplots$x <- 0
      biplots$y <- 0
      biplots$lengths <- sqrt(biplots[,1]^2 + biplots[,2]^2) # Biplot arrow length
      #if(colnames(biplots)[2] == "MDS1"){ # Glitch when only one variable identified? (Need to reverse angle calculation?)
      #  biplots$angles <- atan2(biplots[,2], biplots[,1]) # Biplot arrow angle (radians)
      #}else{
      biplots$angles <- atan2(biplots[,2], biplots[,1]) # Biplot arrow angle (radians)
      #}
      #species <- as.data.frame(scores(rda, display = "sp"))  
      #constraints <- as.data.frame(scores(rda, display = "lc"))

      p <- ggplot(data = sites) + geom_point(aes(x = sites[,1], y = sites[,2], colour = Elevation), shape = 1, size = 2, alpha = 0.6) +
                  geom_text(data = sites, aes(x = sites[,1], y = sites[,2], label = shortIDs, 
                                              colour = Elevation), alpha = 0.6, vjust = 2.5, size = 2) +
                  #geom_polygon(data = sites, aes(x = CAP1, y = CAP2, mapping = Plot, colour = Elevation), alpha = 0.25) +
                  #geom_segment(data = biplots, aes(x = 0, xend = CAP1, y = 0, yend = CAP2),
                  #              arrow = arrow(length = unit(0.25, "cm")), colour="blue") +
                  geom_spoke(data = biplots, aes(x, y, angle = angles, 
                                 radius = lengths*2), arrow = arrow(length = unit(0.25, "cm")), colour = "blue") +
                  geom_text(data = biplots, aes(x = biplots[,1]*2.2, y = biplots[,2]*2.2, 
                                                label = rownames(biplots)), size = 4, colour = "blue", alpha = 0.75) +
                  xlab("") + ylab("") +
                  scale_colour_gradientn(colours = c("blue", "orange")) +
                  theme(panel.grid = element_blank(), plot.title = element_text(size = 8), 
                  plot.margin = unit(c(0,0.25,0.25,0), "cm")) + labs(colour="Elevation (m)") +
                  #coord_fixed() +
                  ggtitle(paste0(letters[n], ". ", names(rda.pick)[[n]], ", Jaccard distance RDA"))
      legend <- get_legend(p)  # Get legend
      
      p <- p + theme(legend.position = "none")
      plot.list[[n]] <- ggplotGrob(p)
      names(plot.list)[[n]] <- names(rda.pick)[[n]]
      n <- n + 1
  }
  return(list(plot.list, legend))
}
        
### Print plots, with size and layout dependent on number of plots ###
output.plots <- function(plots.legend, file.label){
  plot.list <- plots.legend[[1]]
  legend <- plots.legend[[2]]
  ### Split plot.list for printing if length > 6 (too many plots for one page) ###
  if(length(plot.list) == 1){
    pdf(file = paste0("LBI_", file.label, "_min2_jaccard_RDA_", pick, "_plots_revised.pdf"),
        width = 11/2.54, height = 9/2.54, useDingbats = FALSE)
    print (grid.arrange(plot.list[[1]], legend, nrow = 1, ncol=2, widths=c(1, 0.2)))
  }else if(length(plot.list) == 2){
    pdf(file = paste0("LBI_", file.label, "_min2_jaccard_RDA_", pick, "_plots_revised.pdf"),
        width = 20/2.54, height = 9/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=2, nrow=1))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  }else if(length(plot.list) > 2 && length(plot.list) < 5){
    pdf(file = paste0("LBI_", file.label, "_min2_jaccard_RDA_", pick, "_plots_revised.pdf"),
        width = 20/2.54, height = 18/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=2, nrow=2))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  }else if(length(plot.list) > 4 && length(plot.list) < 7){
    pdf(file = paste0("LBI_", file.label, "_min2_jaccard_RDA_", pick, "_plots_revised.pdf"),
        width = 20/2.54, height = 27/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=2, nrow=3))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  }else if(length(plot.list) > 6 && length(plot.list) < 9){
    pdf(file = paste0("LBI_", file.label, "_min2_jaccard_RDA_", pick, "_plots_revised.pdf"),
        width = 20/2.54, height = 36/2.54, useDingbats = FALSE)
    args.list <- c(plot.list, list(ncol=2, nrow=4))
    print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.1)))
  }
  dev.off()
}

###############################################################################
#metric <- "Jaccard"
#metric <- "1-beta1"

genes <- c("16S", "18S")
genes <- c("16S", "18S", "26S", "ITS", "FolCO1", "ShCO1")
taxa <- c("proks-euks")
file.label <- "16S_18S_proks-euks"

df.list <- get.datasets(genes, taxa)
rda.lists <- do.rda(df.list, file.label)
pick <- "Backward"
rda.pick <- rda.lists[[ pick ]]
plots.legend <- get.rda.plots(rda.pick, file.label)
output.plots(plots.legend, file.label)

###############################################################################
# Below: unneeded stuff (old)
OTUtable1 <- t(OTUtable)
OTUtable1 <- OTUtable1[, colSums(OTUtable1) > 0] # Remove empty rows

# Ensure OTUtable and  sample data rownames match, otherwise analysis doesn't work
OTUtable1 <- OTUtable1[order(rownames(OTUtable1)), ]
env_pilot <- env_pilot[order(rownames(env_pilot)), ]

# Constrained ordination ------------------------------------------------------

# Distance-based redundancy analysis, using capscale

# DB-RDA, empty model
rda_0 <- capscale(OTUtable1 ~ 1, env_pilot, distance = "jaccard")

# DB-RDA, maximal model (bad idea - only use for auto model building)
rda_1 <- capscale(OTUtable1 ~ ., env_pilot, distance = "jaccard")
head(summary(rda_1))
anova(rda_1, by = "terms")
plot(rda_1, display = c("wa", "bp"))
# sp = species scores, wa = site scores, bp = biplot arrows, lc = linear constraints 
plot(rda_1, display = c("wa", "bp")) # Note correlation of biplot arrows

# Variance inflation factor - indicates highly correlated variables
vif.cca(rda_1)

# Test each variable individually
capscale(OTUtable1 ~ pH_pilot, env_pilot, distance = "jaccard")
capscale(OTUtable1 ~ Elevation, env_pilot, distance = "jaccard")
capscale(OTUtable1 ~ Mean.Temp, env_pilot, distance = "jaccard")
capscale(OTUtable1 ~ tot_org_carb, env_pilot, distance = "jaccard")
capscale(OTUtable1 ~ tot_nitro, env_pilot, distance = "jaccard")
capscale(OTUtable1 ~ water_content_soil, env_pilot, distance = "jaccard")
capscale(OTUtable1 ~ slope_gradient, env_pilot, distance = "jaccard")
capscale(OTUtable1 ~ Northness, env_pilot, distance = "jaccard")
capscale(OTUtable1 ~ Eastness, env_pilot, distance = "jaccard")

add1(rda_0, scope=formula(rda_1), test="perm")

# Build model automatically (do forward and backward methods give similar results?)
rda_f <- ordistep(rda_0, scope = formula(rda_1), direction = "forward")
rda_b <- ordistep(rda_1, scope = formula(rda_0), direction = "backward")
head(summary(rda_f))
head(summary(rda_b))
anova(rda_f, by = "terms")
anova(rda_b, by = "terms")
plot(rda_f, display = c("wa", "bp"))
plot(rda_b, display = c("wa", "bp"))

# Build model after stepwise removal of collinear variables (vif > 10; requires vif_function.R) 
env_reduced <- vif_func(in_frame = env_pilot)
env_reduced # Remaining variables

# Build model from reduced variable set
# (Unsure how to pass env_reduced variables to capscale formula; paste() doesn't work...)
rda_reduced <- capscale(OTUtable1 ~ pH_pilot + Mean.Temp + tot_nitro + slope_gradient + 
                           Northness+Eastness, env_pilot, distance = "jaccard")
head(summary(rda_reduced))
anova(rda_reduced, by = "terms")
plot(rda_reduced, display = c("wa", "bp"))

# Build model automatically from reduced variable set
rda_reduced_f <- ordistep(rda_0, scope = formula(rda_reduced), direction = "forward")
rda_reduced_b <- ordistep(rda_reduced, scope = formula(rda_0), direction = "backward")
head(summary(rda_reduced_f))
head(summary(rda_reduced_b))
anova(rda_reduced_f, by = "terms")
anova(rda_reduced_b, by = "terms")
plot(rda_reduced_f, display = c("wa", "bp"))
plot(rda_reduced_b, display = c("wa", "bp"))


# Build a model manually
rda_m <- update(rda_0, . ~ . + Mean.Temp + slope_gradient + Northness + Eastness) 
head(summary(rda_m))
anova(rda_m, by = "terms")


