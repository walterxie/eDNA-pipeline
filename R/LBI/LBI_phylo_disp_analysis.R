library(data.table)
library(ggplot2)
library(gridExtra)
library(plyr)
library(reshape2)
library(scales)
library(vegan)
library(pez)
library(picante)
#library(doParallel)
#registerDoParallel(cores = 4)

theme_set(theme_bw(base_size = 8))

#setwd("C:/Documents and Settings/Andrew/Desktop/LBI_miseq_analysis_stuff")
#setwd("C:/Users/Andrew/Desktop/LBI_U8_OTUs_analyses/")
setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/")
#setwd("J:/PhD/PhD Research/NZGL01401_analysis/")

### Load environmnental data ###
envdata.bysubplot <- read.table("LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_subplot.txt", 
                                sep="\t", header=T, row.names=1)
envdata.byplot <- read.table("LBI_U8_OTUs_analyses/LBI_environmental_data/LBI_all_env_data_by_plot.txt", 
                             sep = "\t", header = TRUE, row.names = 1)
#envdata.bysubplot$plot <- rownames(envdata.bysubplot)
envdata.byplot$plot <- rownames(envdata.byplot)

### Path for phylogenetic trees ###
tree_path <- "LBI_U8_OTUs_analyses/LBI_U8_RAxML_trees_v2/FFT-NS2_min2_RAxML_trees/"

### Get correlation statistics ###
get.corr <- function(mpd.mntd, envdata, e){
  cor.list <- list()
  i <- 1
  for(v in unique(mpd.mntd$variable)){
    for(g in unique(mpd.mntd$gene)){
      for(tax in unique(mpd.mntd$taxon)){
        x <- mpd.mntd[(mpd.mntd$variable==v & mpd.mntd$gene==g & mpd.mntd$taxon==tax),]
        if(nrow(x) > 10){ 
          print(paste(v, g, tax, env))
          x$e <- envdata[e][match(as.character(x$plot), envdata$plot),]
          #t <- cor.test(x$obs, x$Elevation, alternative="two.sided", method=c("pearson"))
          t <- cor.test(x$obs, x$e, alternative="two.sided", method=c("pearson"))
          res <- list("variable" = v, "env.var" = e, "gene" = g, "taxon" = tax, "t.stat" = t$statistic, 
                      "pval" = t$p.value, "corr" = t$estimate[1], "df" = t$parameter, 
                      "95_lower" = t$conf.int[1], "95_upper" = t$conf.int[2], "alt" = t$alternative)
          cor.list[[i]] <- res
          i <- i + 1
        }
      }
    }
  }
  return(cor.list)
}


### Loop calculates phylogenetic dispersion (etc.) metrics for all genes and trees ###
genes <- c("16S","18S","26S","ITS","ShCO1","FolCO1")
#genes <- c("26S","ITS","ShCO1","FolCO1")
#g <- "18S"
all.data <- data.frame()
for(g in genes){
  disp.data <- data.frame()
  trees <- Sys.glob(paste0(tree_path,"RAxML_bestTree.", g,"*_OTUs_MAFFT_rxGTR"))
  for(t in trees){
    taxon <- sapply(strsplit(t, "bestTree.*min2_|_OTUs_MAFFT_rxGTR"), "[[", 2)
    print(paste(g, taxon, "..."))
    #t <- read.tree(paste0(tree_path, "RAxML_bestTree.", g, "_LBI_contigs_min2_metazoans_OTUs_MAFFT_rxGTR"))
    tr <- read.tree(t)
    is.rooted(tr)
    f <- Sys.glob(paste0("LBI_U8_OTUs_analyses/OTUtables/", g, "*otutable_min2_by_plot.txt"))
    df <- read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    df.t <- t(df)

    # Make comm object (Warning about dropping columns expected, as OTU tables contain everything but trees don't)
    c <- comparative.comm(tr, df.t, env = envdata.byplot)
    
    #pdiss <- pez.dissimilarity(c) # Phylogenetic distance metrics, returns distance objects 
    pdisp <- pez.dispersion(c, null.model = "independentswap", abundance = FALSE) # Dispersion metrics (slow)
    #pe <- pez.evenness(c) # Phylogenetic structure metrics, includes abundances
    #ps <- pez.shape(c) # Phylogenetic structure metrics, ignores abundances
    pdisp$Plot <- rownames(pdisp)
    pdisp$gene <- g
    pdisp$taxon <- taxon
    disp.data <- rbind(disp.data, pdisp)
  }
  write.table(disp.data, file = paste0(g, "_min2_phylo_dispersion_data.txt"), sep = "\t", quote = FALSE, col.names = NA)
  all.data <- rbind(all.data, disp.data)
}

write.table(all.data, file = "All_genes_min2_phylo_dispersion_data.txt", sep = "\t", quote = FALSE, col.names = NA)

all.data <- read.table(paste0("LBI_U8_OTUs_analyses/Phylo_dispersion_analysis/All_genes_min2_phylo_dispersion_data.txt"), 
                       sep = "\t", header = TRUE, row.names = 1)
#colnames(all.data) <- gsub("plot", "Plot", colnames(all.data))
all.data$gene <- gsub("FolCO1", "COI-650", all.data$gene)
all.data$gene <- gsub("ShCO1", "COI-300", all.data$gene)

### Organise mpd and mntd data for plotting ###
mpd <- all.data[, c("ses.mpd.mpd.obs.z","ses.mpd.mpd.obs.p","gene","taxon","plot")]  
mpd$variable <- "ses.mpd"
colnames(mpd)[1] <- "obs"
colnames(mpd)[2] <- "p"
mpd$p.adj <- p.adjust(mpd$p, method = "BH")

mntd <- all.data[, c("ses.mntd.mntd.obs.z","ses.mntd.mntd.obs.p","gene","taxon","plot")]
mntd$variable <- "ses.mntd"
colnames(mntd)[1] <- "obs"
colnames(mntd)[2] <- "p"
mntd$p.adj <- p.adjust(mntd$p, method = "BH")

mpd.mntd <- rbind(mpd, mntd)

### Organise for different symbols to indicate significance values ###
mpd.mntd$sig <- mpd.mntd$p.adj 
mpd.mntd$sig[mpd.mntd$p.adj  <= 0.05] <- "S"
mpd.mntd$sig[mpd.mntd$p.adj  > 0.05] <- "NS"
mpd.mntd$obs.1 <- -1*mpd.mntd$obs

mpd.mntd$gene <- factor(mpd.mntd$gene, levels = c("16S","18S","26S","ITS","COI-300","COI-650"), ordered = TRUE)  
mpd.mntd$taxon <- factor(mpd.mntd$taxon, levels = c("bacteria","protists","fungi","metazoans","plants"), ordered = TRUE)  
mpd.mntd$variable <- factor(mpd.mntd$variable, levels = c("ses.mpd","ses.mntd"), ordered = TRUE)
mpd.mntd$variable2 <- mpd.mntd$variable
mpd.mntd$variable2 <- gsub("ses.mpd","Nearest relative index", mpd.mntd$variable2)
mpd.mntd$variable2 <- gsub("ses.mntd","Net taxon index", mpd.mntd$variable2)
mpd.mntd$variable2 <- factor(mpd.mntd$variable2, levels = c("Nearest relative index","Net taxon index"), ordered = TRUE)

#mpd.mntd$Elevation <- envdata.byplot$Elevation[match(mpd.mntd$plot, rownames(envdata.byplot))]
#mpd.mntd$pH <- envdata.byplot$pH[match(mpd.mntd$plot, rownames(envdata.byplot))]
envdata <- envdata.byplot[, c(1:4,9:16,23)]
### Log transform skewed chemical variables ###
envdata$log.EC <- log(envdata$EC)
envdata$log.Organic.C <- log(envdata$Organic.C) 
envdata$log.Total.N <- log(envdata$Total.N)
envdata$log.NO3.N <- log(envdata$NO3.N)
envdata$log.NH4.N <- log(envdata$NH4.N)
envdata$log.Olsen.P <- log(envdata$Olsen.P)

#mpd.mntd <- merge(mpd.mntd, envdata) # Combine with all env data?

### Get correlation statistics ###
#env <- "Elevation"
#cor.data <- rbindlist(get.corr(mpd.mntd, env))
#write.table(cor.data, file = paste0("All_genes_min2_phylo_disp_corr_data.txt"), sep = "\t", quote = FALSE)
#cor.data <- read.table(paste0("LBI_U8_OTUs_analyses/Phylo_dispersion_analysis/All_genes_min2_phylo_dispersion_data.txt"), 
#                       sep = "\t", header = TRUE, row.names = 1)

### Get correlations for different environmental variables ###
all.corr.data <- data.frame()
env.vars <- c("Elevation","pH","log.EC","log.Organic.C","log.Total.N","C.N.ratio","log.NO3.N","log.NH4.N","log.Olsen.P")
for(e in env.vars){
  print(e)
  cor.data <- rbindlist(get.corr(mpd.mntd, envdata, e))
  cor.data$p.adj <- p.adjust(cor.data$pval, method = "BH")
  #write.table(cor.data, file = paste0("All_genes_min2_phylo_disp_corr_data_", env, ".txt"), sep = "\t", quote = FALSE)
  all.corr.data<- rbind(all.corr.data, cor.data)
}

### Adjust overall p-values (?) and write to file ### 
mpd.corr <- all.corr.data[all.corr.data$variable == "ses.mpd", ]
mntd.corr <- all.corr.data[all.corr.data$variable == "ses.mntd", ]
mpd.corr$ovr.p.adj <- p.adjust(mpd.corr$pval, method = "BH")
mntd.corr$ovr.p.adj <- p.adjust(mntd.corr$pval, method = "BH")
write.table(mpd.corr, file = "All_genes_mpd_all_correlation_data.txt", sep = "\t", quote = FALSE, col.names = NA)
write.table(mntd.corr, file = "All_genes_mntd_all_correlation_data.txt", sep = "\t", quote = FALSE, col.names = NA)
all.corr.data <- rbind(mpd.corr, mntd.corr)

### Drop some datasets for plotting ###
mpd.mntd <- mpd.mntd[!(mpd.mntd$gene == "ITS" & mpd.mntd$taxon == "protists"), ]
mpd.mntd <- mpd.mntd[!(mpd.mntd$gene == "COI-650" & mpd.mntd$taxon == "fungi"), ]
#pd.mntd.2 <- mpd.mntd.2[!(mpd.mntd.2$gene == "18S" & mpd.mntd.2$taxon == "plants"), ]
#mpd.mntd.2 <- mpd.mntd.2[!(mpd.mntd.2$gene == "26S" & mpd.mntd.2$taxon == "plants"), ]

for(e in env.vars){
### Make a nice figure ###
  mpd.mntd$e <- envdata[e][match(as.character(mpd.mntd$plot), envdata$plot),]
  p <- ggplot(mpd.mntd, aes(x = e, y = obs.1)) + 
    geom_point(data = mpd.mntd, aes(colour = taxon, shape = sig)) + 
    scale_shape_manual(values = c(1,19), guide=FALSE) +
    #geom_smooth(data = mpd.mntd, method = "lm", mapping = aes(x = Elevation, y = -1*obs, colour = taxon)) +
    ylab("") +
    geom_hline(yintercept = 0, colour = "grey") +
    facet_grid(variable ~ gene, scales = "free") +
    scale_color_brewer(palette = "Set1") +
    theme(panel.grid = element_blank(), strip.background = element_blank())
  
  ### Make a nice figure with correlations added ###
  cor.05 <- all.corr.data[all.corr.data$env.var == e & all.corr.data$ovr.p.adj < 0.05,] # get "significant" correlations only
  mpd.mntd$cor.05 <- cor.05$ovr.p.adj[match(interaction(mpd.mntd$gene, mpd.mntd$taxon, mpd.mntd$variable), 
                           interaction(cor.05$gene, cor.05$taxon, cor.05$variable))]
  
  #p + geom_line(data = mpd.mntd, # Include all correlations
  #               mapping = aes(x = Elevation, y = obs.1, colour = taxon), 
  #               alpha = 0.75, size = 0.7, stat = "smooth", method = "lm") 
  
  p <- p + geom_line(data = mpd.mntd[!(is.na(mpd.mntd$cor.05)), ], # Includes significant correlations only
                  mapping = aes(x = e, y = obs.1, colour = taxon), 
                  alpha = 0.75, size = 0.7, stat = "smooth", method = "lm")
  
  ggsave(p, file = paste0("All_genes_min2_phylo_disp_corr_", e, "_padj.pdf"), 
         width = 297, height = 210, units = "mm", useDingbats = FALSE)
}

###############################################################################

all.corr.data$gene_taxon <- paste(all.corr.data$gene, all.corr.data$taxon)
all.corr.data <- all.corr.data[all.corr.data$gene_taxon != "ITS protists", ]
all.corr.data <- all.corr.data[all.corr.data$gene_taxon != "COI-650 fungi", ]
all.corr.data$gene_taxon <- factor(all.corr.data$gene_taxon, 
                                   levels = rev(c("16S bacteria","18S protists","18S fungi","18S metazoans","18S plants",
                                                  "26S protists","26S fungi","26S metazoans","26S plants","ITS fungi",
                                                  "COI-300 protists","COI-300 fungi","COI-300 metazoans","COI-300 plants",
                                                  "COI-650 bacteria","COI-650 protists","COI-650 fungi","COI-650 metazoans",
                                                  "COI-650 plants")), ordered = TRUE)
all.corr.data$env.var <- factor(all.corr.data$env.var, 
                                levels = c("Elevation","pH","log.EC","log.Organic.C","log.Total.N","log.NO3.N",
                                           "log.NH4.N","C.N.ratio","log.Olsen.P"), ordered = TRUE)  
  
#all.corr.data$corr_sig <- paste0(round(all.corr.data$corr, 2), "\n(", all.corr.data$sig, ")")
#all.corr.data$corr_sig <- gsub("0.001", "<0.001", all.corr.data$corr_sig)
all.corr.data$sig_ind <- all.corr.data$ovr.p.adj
all.corr.data$sig_ind <- ifelse(all.corr.data$ovr.p.adj <= 0.001, "***", all.corr.data$ovr.p.adj)
all.corr.data$sig_ind <- ifelse(all.corr.data$sig_ind > 0.001 & all.corr.data$sig_ind <= 0.01, "**", all.corr.data$sig_ind)
all.corr.data$sig_ind <- ifelse(all.corr.data$sig_ind > 0.01 & all.corr.data$sig_ind <= 0.05, "*", all.corr.data$sig_ind)
all.corr.data$sig_ind <- ifelse(all.corr.data$sig_ind > 0.05, "", all.corr.data$sig_ind)  
  
ggplot(all.corr.data, aes(x = env.var, y = gene_taxon, fill = corr)) + geom_tile() +
    facet_grid(. ~ variable) +
    #geom_text(data = d, aes(x = t1, y = t2, label = d$corr_sig), size = 1.8) +
    geom_text(data = all.corr.data, aes(x = env.var, y = gene_taxon, label = sig_ind), size = 1.8) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7), 
          axis.text.y = element_text(size = 7), panel.grid = element_blank(), 
          legend.title = element_text(size = 7), legend.key.width = unit(0.5, "cm"),
          strip.background = element_rect(fill=NA, colour = NA)) +
    scale_fill_gradient2(high = "#f46d43", mid = "#ffffbf", low = "#3288bd",
                         midpoint = 0, limit = c(-1,1), name = "Correlation") +
    #  scale_fill_gradient2(high = "#f1a340", mid = "#f7f7f7", low = "#998ec3",
    #                       midpoint = 0.5, limit = c(0,1), name = "Correlation") +
    #  scale_fill_gradient(high = "#f46d43", low = "#3288bd", limit = c(0,1), name = "Correlation") +
    xlab("") + ylab("")
  
