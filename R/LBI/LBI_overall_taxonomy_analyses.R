library(data.table)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)
library(vegetarian)

theme_set(theme_bw(base_size=8))

# Colour palette
library(RColorBrewer)
colors <- brewer.pal(8, "Spectral")
#colors <- brewer.pal(4, "Paired")
#colors <- brewer.pal(4, "Set3")
pal <- colorRampPalette(colors) 

#tax_ref <- read.table("C:/Documents and Settings/Andrew/Desktop/Database stuff/New_taxonomy_from_PLOSONE_2015_fixed.txt", 
#                      header = TRUE, sep = "\t", quote = "", comment.char = "")
tax_ref <- read.table("H:/My Documents/PhD Research PFR folder/Database stuff/New_taxonomy_from_PLOSONE_2015_cleaned_up.txt", 
                      header = TRUE, sep = "\t", quote = "", comment.char = "")

#setwd("C:/Documents and Settings/Andrew/Desktop/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/")
setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/")

### Load environmnental data ###
envdata.bysubplot <- read.table("LBI_environmental_data/LBI_all_env_data_by_subplot.txt", sep="\t", header=T, row.names=1)
envdata.byplot <- read.table("LBI_environmental_data/LBI_all_env_data_by_plot.txt", 
                             sep = "\t", header = TRUE, row.names = 1)
envdata.bysubplot$Plot <- rownames(envdata.bysubplot)
envdata.byplot$Plot <- rownames(envdata.byplot)

### Load OTU datasets ###
get.datasets <- function(genes){
  df.list <- list()
  n <- 1
  for(g in genes){
    #print(g)
    f <- Sys.glob(paste0("OTUtables/", g, "*otutable.txt"))
    #f <- Sys.glob(paste0("OTUtables/", g, "*otutable_by_plot.txt"))
    df <- read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    t <- Sys.glob(paste0("Taxa_tables/", g, "*s_nt_paths_new_taxonomy_table.txt"))
    taxonomy <- read.table(t, sep="\t", header=TRUE, row.names=1, quote = "", comment.char = "")
    ### Remove assorted quirks in taxonomy! ###
    pattern.x <- "(\\s\\[|\\()(\\=|\\.|\\,|\\s|\\w|\\?)*(\\]|\\))" 
    taxonomy$Kingdom <- gsub(pattern.x, "", taxonomy$Kingdom, perl = TRUE)
    taxonomy$Phylum <- gsub(pattern.x, "", taxonomy$Phylum, perl = TRUE)
    taxonomy$Class <- gsub(pattern.x, "", taxonomy$Class, perl = TRUE)
    taxonomy$Order <- gsub(pattern.x, "", taxonomy$Order, perl = TRUE)
    df$Superkingdom <- taxonomy$Superkingdom[match(rownames(df), rownames(taxonomy))]
    df$Kingdom <- taxonomy$Kingdom[match(rownames(df), rownames(taxonomy))]
    df$Phylum <- taxonomy$Phylum[match(rownames(df), rownames(taxonomy))]
    df$Class <- taxonomy$Class[match(rownames(df), rownames(taxonomy))]
    df$Order <- taxonomy$Order[match(rownames(df), rownames(taxonomy))]
    if(g == "ShCO1"){ # Fix COI names
      g <- "COI-300"
    }
    if(g == "FolCO1"){
      g <- "COI-650"
    }
    df.list[[n]] <- df
    label <- paste(g)
    print(label)
    print(nrow(df))
    names(df.list)[[n]] <- label
    n <- n + 1
  }
  return(df.list)
}
    
### Get OTU and sequence counts by taxonomic group ###
get.counts.sums <- function(df.list, tax.rank){
  counts.sums.list <- list()
  n <- 1
  for(df in df.list){
    OTU.sums <- cbind(df[,57:61], rowSums(df[,1:56]))
    #OTU.sums <- cbind(df[,29:33], rowSums(df[,1:28]))
    colnames(OTU.sums)[6] <- "OTU.sum"

    dt <- data.table(OTU.sums)
    counts.1 <- dt[OTU.sum > 0, .("count.1" = .N), by = tax.rank]
    counts.2 <- dt[OTU.sum > 1, .("count.2" = .N), by = tax.rank]
    sums.1 <- dt[OTU.sum > 0, .("sum.1" = sum(OTU.sum)), by = tax.rank]
    sums.2 <- dt[OTU.sum > 1, .("sum.2" = sum(OTU.sum)), by = tax.rank]
    
    counts.sums <- merge(counts.1, sums.1, by = tax.rank)
    counts.sums.2 <- merge(counts.2, sums.2, by = tax.rank)
    counts.sums <- merge(counts.sums, counts.sums.2, by = tax.rank, all = TRUE)
    counts.sums[is.na(counts.sums)] <- 0
    counts.sums$gene <- paste(names(df.list)[[n]])

    label <- paste(names(df.list)[[n]], tax.rank)
    print(label)
    print(nrow(counts.sums))
    counts.sums.list[[n]] <- counts.sums
    names(counts.sums.list)[[n]] <- label
    n <- n + 1
  }
  return(rbindlist(counts.sums.list))
}

### Get biodiversity stats by taxonomic group ###
get.div.by.taxa <- function(df.list, tax.rank, m){
  all.div <- list()
  n <- 1
  for(df in df.list){
    div.list <-  list()
    i <- 1
    
    df.1 <- df[!is.na(df$Phylum), ] # Min 1 stats
    if(m == "min2"){
      df.1 <- df.1[rowSums(df.1[,1:56]) > 1, ] # Min 2 stats 
      #df.1 <- df.1[rowSums(df.1[,1:28]) > 1, ]
    }
    for(lev in unique(df.1[, tax.rank])){
        print(lev)
        tax.subset.1 <- df.1[grep(lev, df.1[, tax.rank]), 1:56]
        #tax.subset.1 <- df.1[grep(lev, df.1[, tax.rank]), 1:28]
        print(dim(tax.subset.1))
        if(nrow(tax.subset.1) > 0){
          a0 <- d(t(tax.subset.1), lev = "alpha", q = 0)
          a1 <- d(t(tax.subset.1), lev = "alpha", q = 1)
          a2 <- d(t(tax.subset.1), lev = "alpha", q = 2)
          b0 <- d(t(tax.subset.1), lev = "beta", q = 0)
          b1 <- d(t(tax.subset.1), lev = "beta", q = 1)
          b2 <- d(t(tax.subset.1), lev = "beta", q = 2)
          g0 <- d(t(tax.subset.1), lev = "gamma", q = 0)
          g1 <- d(t(tax.subset.1), lev = "gamma", q = 1)
          g2 <- d(t(tax.subset.1), lev = "gamma", q = 2)
          div <- list("min" = m, "tax.rank" = lev, "gene" = names(df.list)[[n]], 
                        "alpha0" = a0, "alpha1" = a1, "alpha2" = a2,
                        "beta0" = b0, "beta1" = b1, "beta2" = b2,
                        "gamma0" = g0, "gamma1" = g1, "gamma2" = g2)
          div.list[[i]] <- div
          i <- i + 1
        }
    }
    all.div[[n]] <- div.list 
    n <- n + 1
  }
  return(rbindlist(all.div))
}  


### Get biodiversity stats by taxonomic group and plot ###
get.div.by.taxa.plot <- function(df.list, tax.rank, m){
  all.div <- list()
  n <- 1
  for(df in df.list){
    div.list <-  list()
    i <- 1
    
    df.1 <- df[!is.na(df$Phylum), ] # Min 1 stats
    if(m == "min2"){
      df.1 <- df.1[rowSums(df.1[, 1:56]) > 1, ] # Min 2 stats 
      #df.1 <- df.1[rowSums(df.1[,1:28]) > 1, ]
    }
    
    plots <- gsub("\\w$", "", colnames(df.1[, 1:56]))
    for(lev in unique(df.1[, tax.rank])){
        print(lev)
        tax.subset <- df.1[grep(lev, df.1[, tax.rank]), ]
        for(pl in unique(plots)){
          print(pl)
          df.2 <- tax.subset[, grep(pl, colnames(tax.subset))]
          print(dim(df.2))
          if(nrow(df.2) > 0){
            a0 <- d(t(df.2), lev = "alpha", q = 0)
            a1 <- d(t(df.2), lev = "alpha", q = 1)
            a2 <- d(t(df.2), lev = "alpha", q = 2)
            b0 <- d(t(df.2), lev = "beta", q = 0)
            b1 <- d(t(df.2), lev = "beta", q = 1)
            b2 <- d(t(df.2), lev = "beta", q = 2)
            g0 <- d(t(df.2), lev = "gamma", q = 0)
            g1 <- d(t(df.2), lev = "gamma", q = 1)
            g2 <- d(t(df.2), lev = "gamma", q = 2)
            div <- list("min" = m, "tax.rank" = lev, "gene" = names(df.list)[[n]], "plot" = gsub("-", "", pl),
                          "alpha0" = a0, "alpha1" = a1, "alpha2" = a2,
                          "beta0" = b0, "beta1" = b1, "beta2" = b2,
                          "gamma0" = g0, "gamma1" = g1, "gamma2" = g2)
            div.list[[i]] <- div
            i <- i + 1
        }
      }
    }
    all.div[[n]] <- rbindlist(div.list) 
    n <- n + 1
  }
  return(rbindlist(all.div))
} 


### Fix x-axis number format for plot ###
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

### Make a plot of OTU and sequence counts by taxonomic group ###
make.taxonomy.plot <- function(all.counts.sums){
  z <- all.counts.sums[all.counts.sums$Phylum != 0, ]
  z$Kingdom <- tax_ref$Kingdom[match(z$Phylum, tax_ref$Phylum)]
  z <- melt(z, id.vars = c("Kingdom","Phylum","gene"))
  # Order factors
  z$Phylum <- factor(z$Phylum, ordered = TRUE, levels = rev(tax_ref$Phylum))
  z$gene <- factor(z$gene, ordered = TRUE, levels = c("16S", "18S", "26S", "ITS", "COI-300", "COI-650"))
  z$Kingdom <- gsub("root|No hits|Not assigned|cellular organisms", "Unknown", z$Kingdom)
  z$Kingdom <- factor(z$Kingdom, ordered = TRUE, levels = c("ARCHAEA","BACTERIA","EUKARYOTA","PROTOZOA",
                                                            "CHROMISTA","FUNGI","PLANTAE","ANIMALIA","Unknown"))
  p <- ggplot(na.omit(z)) + 
    geom_point(data = z[variable == "count.1"], aes(x = Phylum, y = value, colour = Kingdom), shape = 1, size = 2) +
    geom_point(data = z[variable == "sum.1"], aes(x = Phylum, y = value, colour = Kingdom), shape = 5, size = 2, show_guide = F) +
    geom_point(data = z[variable == "count.2"], aes(x = Phylum, y = value, colour = Kingdom), shape = 16, size = 2) +
    #geom_point(data = z[variable == "sums.2"], aes(x = Phylum, y = value, colour = group1), shape = 17, size = 2) +
    geom_line(data = na.omit(z), aes(x = Phylum, y = value, group=interaction(Phylum, gene), 
                                     colour = Kingdom), size = 0.5, alpha = 0.5) + 
    facet_grid( ~ gene) + coord_flip() + guides(fill = guide_legend(reverse = FALSE)) +
    ylab("Number of sequences or OTUs") + xlab("Phylum (or higher-level taxon)") + 
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), label=scientific_10) +
    #scale_color_manual(values = pal(length(unique(z$Kingdom)))) +
    theme(strip.background = element_blank(), plot.title = element_text(size = 9))
  return(p)
}

###############################################################################
genes <- c("16S", "18S", "26S", "ITS", "FolCO1", "ShCO1")
tax.rank <- "Phylum" # Rank to summarise data by
tax.rank <- "Kingdom"
m <- "min2"

df.list <- get.datasets(genes)

### Make some adjustments to taxonomy? ###
i <- 1
for(df in df.list){
  print(dim(df))
  df$Kingdom <- gsub("CHROMISTA|PROTOZOA", "PROTISTS", df$Kingdom)
  df.list[[i]] <- df
  i <- i + 1 
}

all.counts.sums <- get.counts.sums(df.list, tax.rank)
div.by.taxa <- get.div.by.taxa(df.list, tax.rank, m)
div.by.taxa.plot <- get.div.by.taxa.plot(df.list, tax.rank, m)

all.counts.sums$Kingdom <- tax_ref$Kingdom[match(all.counts.sums$Phylum, tax_ref$Phylum)]
all.counts.sums$Kingdom <- gsub("root|No hits|Not assigned|cellular organisms", "Unknown", all.counts.sums$Kingdom)
write.table(all.counts.sums, file = paste0("Overall_counts_sums_by_", tax.rank, ".txt"), 
            sep = "\t", quote = FALSE, col.names = NA)

all.div$Kingdom <- tax_ref$Kingdom[match(all.div$tax.rank, tax_ref$Phylum)]
all.div$Kingdom <- gsub("root|No hits|Not assigned|cellular organisms", "Unknown", all.div$Kingdom)
write.table(div.by.taxa, file = paste0("Overall_alpha_div_stats_by_plot_by_", tax.rank, ".txt"), 
            sep = "\t", quote = FALSE, col.names = NA)

p1 <- make.taxonomy.plot(all.counts.sums)
ggsave(p1, file = "LBI_overall_taxonomy_OTUs_reads_by_Phylum.pdf", width = 260, height = 200, units = "mm")


write.table(div.by.taxa.plot, file = "LBI_biodiversity_stats_by_kingdom-plot.txt", sep = "\t", quote = FALSE, col.names = NA)

### Drop various gene/taxon combinations for table ###
d2 <- div.by.taxa.plot
d2 <- d2[!(grepl("root|cellular organisms|Not assigned|No hits", d2$tax.rank)), ]
d2 <- d2[!(d2$gene == "16S" & d2$tax.rank != "BACTERIA"), ]
d2 <- d2[!(grepl("18S|26S|ITS|COI-300|COI-650", d2$gene) & d2$tax.rank == "ARCHAEA"), ]
d2 <- d2[!(grepl("18S|26S|ITS|COI-300", d2$gene) & d2$tax.rank == "BACTERIA"), ]
d2 <- d2[!(grepl("ITS", d2$gene) & d2$tax.rank == "PLANTAE"), ]
d2 <- d2[!(grepl("EUKARYOTA", d2$tax.rank)), ]

vars <- c("alpha0","alpha1","alpha2","beta0","beta1","beta2","gamma0","gamma1","gamma2")
for(v in vars){
  d3 <- dcast(d2, plot ~ gene + tax.rank, value.var = v)
  d3$Elevation <- envdata.byplot$Elevation[match(d3$plot, envdata.byplot$Plot)]
  d3 <- d3[order(d3$Elevation), ]
  write.table(d3, file = paste0("LBI_", v, "_biodiversity_stats_by_kingdom_table.txt"), sep = "\t", quote = FALSE, col.names = NA)
}

colnames(d3) <- gsub("_", "\n", colnames(d3))
plot(d3[, 2:23], gap = 0, panel = panel.smooth) 
       