# Figure S1: Correlations between spatial separation and multivariate dissimilarity index (DSI)
# Figure S2: Correlations between elevation differences and multivariate dissimilarity index (DSI)
# Figure S3: Within-plot (red) and between-plot (blue) pairwise sample similarities

# do not have within between, if by.plot=T
# all.dist.subplot <- getDissimVsDistances(save.rdata=TRUE)
# all.dist.plot <- getDissimVsDistances(by.plot=T, save.rdata=TRUE)
# Use "save.rdata=T", if you want to save it to "all.dist.list".
getDissimVsDistances <- function(input.names, by.plot=FALSE, save.rdata=FALSE) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)
  output.names <- getOutputNames(input.names)
  
  metrics <- c("jaccard", "horn.morisita", "bray.curtis", "beta1-1", "wt.unif", "unwt.unif")
  
  # calculate all dist
  cat("\nCalculate all dissimilarity ... \n")
  all.dist.list <- list()
  for (data.id in 1:length(input.names)) {
    cm <- getIdentifiedCM(input.names[data.id], by.plot=by.plot)
    tree <- getPhyloTree(input.names[data.id])
    
    for(metric in metrics){
      dsi.dist <- ComMA::getDissimilarity(cm, tree=tree, method=metric)
      all.dist.list[[output.names[data.id]]][[metric]] <- dsi.dist
    }
  }
  all.dist.list$by.plot <- by.plot
#  if (save.rdata)
#    save(all.dist.list, file = paste("data/all.dist", ifelse(by.plot, "plot", "subplot"), "RData", sep = ".") )
  
  all.dist.list <- appendDistDF(output.names, metrics, all.dist.list, by.plot=by.plot)
  if (save.rdata)
    save(all.dist.list, file = paste("data/all.dist", ifelse(by.plot, "plot", "subplot"), "RData", sep = ".") )
  
  return(all.dist.list) # a list of list
}

# do not have within between, if by.plot=T
appendDistDF <- function(output.names, metrics, all.dist.list=list(), by.plot=FALSE) {
  if (!by.plot) {
    # all within-between pairwise distances
    cat("\nCalculate all within-between pairwise distances ... \n")
    dist.data <- data.frame(check.names = F)
    for (data.id in 1:length(output.names)) {
      for(metric in metrics){
        d <- getWithinBetweenDistDF(all.dist.list[[output.names[data.id]]][[metric]])
        if(nrow(d) > 0){
          d$gene <- output.names[data.id]
          d$metric <- metric
        }
        dist.data <- rbind(dist.data, d)
      }
    }
    all.dist.list$within.between <- dist.data
    
    # within-plot spatial distances
    cat("\nCalculate within-plot spatial distances ... \n")
    dist.data <- data.frame(check.names = F)
    for (data.id in 1:length(output.names)) {
      for(metric in metrics){
        d <- getWithinDistDF(all.dist.list[[output.names[data.id]]][[metric]])
        if(nrow(d) > 0){
          d$gene <- output.names[data.id]
          d$metric <- metric
        }
        dist.data <- rbind(dist.data, d)
      }
    }
    all.dist.list$within <- dist.data
    
    # between-plot spatial distances
    cat("\nCalculate between-plot spatial distances ... \n")
    dist.data <- data.frame(check.names = F)
    for (data.id in 1:length(output.names)) {
      for(metric in metrics){
        d <- getBetweenDistDF(all.dist.list[[output.names[data.id]]][[metric]])
        if(nrow(d) > 0){
          d$gene <- output.names[data.id]
          d$metric <- metric
        }
        dist.data <- rbind(dist.data, d)
      }
    }
    all.dist.list$between <- dist.data
  } else {
    cat("\nDo not have within between measurements, if by.plot=TRUE.\n")
  }
  
  # elevation differences between subplots
  cat("\nCalculate elevation differences between subplots ... \n")
  env <- getEnvData(by.plot=by.plot)
  dist.data <- data.frame(check.names = F)
  for (data.id in 1:length(output.names)) {
    for(metric in metrics){
      d <- getElevationDistDF(all.dist.list[[output.names[data.id]]][[metric]], env)
      if(nrow(d) > 0){
        d$gene <- output.names[data.id]
        d$metric <- metric
      }
      dist.data <- rbind(dist.data, d)
    }
  }
  all.dist.list$elev.diff <- dist.data
  return(all.dist.list)
}

renderMetrics <- function(data, data.levels, method.levels) {
  data$gene <- factor(data$gene, levels = data.levels, ordered = TRUE)
  data$metric <- gsub("jaccard","Jaccard", data$metric)
  data$metric <- gsub("bray.curtis","Bray-Curtis", data$metric)
  data$metric <- gsub("horn.morisita","Morisita-Horn", data$metric)
  data$metric <- gsub("unwt.unif","unweighted Unifrac", data$metric)
  data$metric <- gsub("wt.unif","weighted Unifrac", data$metric)
  data$metric <- factor(data$metric, levels = method.levels, ordered = TRUE)
  return(data)
}

# plotDistanceCorrelation(all.dist.list[["within"]])
# plotDistanceCorrelation(all.dist.list[["elev.diff"]])
plotDistanceCorrelation <- function(dist.data, data.levels = c("16S","18S","26S","ITS","COI-300","COI-650"), 
                                    method.levels = c("Jaccard","Bray-Curtis","Morisita-Horn","beta1-1",
                                                      "unweighted Unifrac","weighted Unifrac"),
                                    x.lab="Spatial distance (m)", y.lab="Distance measure") {
  
  dist.data <- renderMetrics(dist.data, data.levels, method.levels)
  
  ### Distance correlation plots ###
  require(ggplot2)
  theme_set(theme_bw(base_size=8))
  p <- ggplot(dist.data, aes(x = dist, y = value)) + 
    geom_point(shape = 1) + 
    geom_smooth(method = "lm", se = FALSE) + 
    ylab(y.lab) + xlab(x.lab) +
    facet_grid(metric ~ gene, scales = "free") +
    theme(strip.background = element_blank(), panel.grid = element_blank())
  
  ### Add correlation data to plot? ###  
  corlist <- getCorrDist(dist.data)
  require(data.table)
  cor.data <- rbindlist(corlist)
  cor.data$cp1 <- paste0(round(cor.data$corr, 2), " (", round(cor.data$pval, 3), ")")
  cor.data$cp <- gsub("\\(0\\)", "\\(<0.001\\)", cor.data$cp1)
  
  cor.data <- renderMetrics(cor.data, data.levels, method.levels)
  
  #long cut way to find number of facets
  len <- length(levels(dist.data$metric))*length(levels(dist.data$gene))
  vars <- data.frame(expand.grid(levels(dist.data$metric), levels(dist.data$gene)))
  colnames(vars) <- c("metric", "gene")
  # x and y are coordinates relative to x/y axis scales
  dat <- data.frame(x = rep(10, len), y = rep(0.5, len), vars) 
  dat$labs <- cor.data$cp[match(interaction(cor.data$metric, cor.data$gene), interaction(dat$metric, dat$gene))]
  
  p <- p + geom_text(data = dat, aes(x, y, label=labs, group=NULL), size = 2)
  return(p)
}

# plotWithinBetween(all.dist.list[["within.between"]])
plotWithinBetween <- function(dist.data, data.levels = c("16S","18S","26S","ITS","COI-300","COI-650"), 
                              method.levels = c("Jaccard","Bray-Curtis","Morisita-Horn","beta1-1",
                                                "unweighted Unifrac","weighted Unifrac"),
                              wb.levels = c("within", "between"),
                              x.lab="", y.lab="Distance measure") {
  if (! "wb" %in% colnames(dist.data) )
    stop("Require all within-between pairwise distance !\nUse getWithinBetweenDistDF.")
  
  dist.data <- renderMetrics(dist.data, data.levels, method.levels)
  dist.data$wb <- factor(dist.data$wb, levels = wb.levels, ordered = TRUE)
  
  ### Within-between distances boxplot ###
  require(ggplot2)
  theme_set(theme_bw(base_size=8))
  p <- ggplot(dist.data, aes(x = wb, y = value, color = wb)) + 
    geom_boxplot(outlier.shape = 1, outlier.colour = alpha("black", 0.25)) + 
    ylab(y.lab) + xlab(x.lab) +
    facet_grid(metric ~ gene, scales = "free") +
    scale_colour_brewer(palette = "Set1") +
    theme(strip.background = element_blank(), panel.grid = element_blank(),
          legend.position = "none")
  return(p)
}


## Get all within-between pairwise distances ###
getWithinBetweenDistDF <- function(dsi.dist){
  require(reshape2)
  dist.list <- melt(as.matrix(dsi.dist))
  dist.list$p1 <- gsub("-.", "", dist.list$Var1)
  dist.list$p2 <- gsub("-.", "", dist.list$Var2)
  dist.list$s1 <- sapply(strsplit(as.character(dist.list$Var1), "-"), "[[", 2)
  dist.list$s2 <- sapply(strsplit(as.character(dist.list$Var2), "-"), "[[", 2)
  dist.list$wb <- ifelse(dist.list$p1 == dist.list$p2, "within", "between")
  dist.list <- dist.list[dist.list$Var1 != dist.list$Var2, ] # Drop same sample distances
  return(dist.list)
}

### Get within-plot spatial distances ###
getWithinDistDF <- function(dsi.dist){
  require(reshape2)
  dist.list <- melt(as.matrix(dsi.dist))
  dist.list$p1 <- gsub("-.", "", dist.list$Var1)
  dist.list$p2 <- gsub("-.", "", dist.list$Var2)
  dist.list$s1 <- sapply(strsplit(as.character(dist.list$Var1), "-"), "[[", 2)
  dist.list$s2 <- sapply(strsplit(as.character(dist.list$Var2), "-"), "[[", 2)
  dist.list <- dist.list[dist.list$p1 == dist.list$p2, ] # Drop between-plot distances
  dist.list <- dist.list[dist.list$s1 != dist.list$s2, ] # Drop same sample distances
  dist.list$spair <- subplotSort(paste0(dist.list$s1, dist.list$s2))
  
  subplot.dists  <- getWithinPlotDistance(verbose=FALSE)
  dist.list$dist <- subplot.dists$dist[match(dist.list$spair, subplot.dists$pair)]
  return(dist.list)
}

### Get between-plot spatial distances ###
getBetweenDistDF <- function(dsi.dist){
  require(reshape2)
  dist.list <- melt(as.matrix(dsi.dist))
  dist.list$p1 <- gsub("-.", "", dist.list$Var1)
  dist.list$p2 <- gsub("-.", "", dist.list$Var2)
  dist.list <- dist.list[dist.list$p1 != dist.list$p2, ] # Drop within-plot distances
  dist.list$pair <- plotSort(paste(dist.list$p1, dist.list$p2)) 
  
  plot.dists <- getBetweenPlotDistance(verbose=FALSE)
  dist.list$dist <- plot.dists$dist[match(dist.list$pair, plot.dists$pair)]
  return(dist.list)
}

### Get elevation differences between plots/subplots ###
# presume rownames(env) are plots/subplots
getElevationDistDF <- function(dsi.dist, env, elev.name="Elevation"){
  if (! elev.name %in% colnames(env) )
    stop("Cannot find elevation column from env file ! ")
  require(reshape2)
  dist.list <- melt(as.matrix(dsi.dist))
  #dist.list$e1 <- env$Elevation[match(dist.list$Var1, env$Plot)]
  #dist.list$e2 <- env$Elevation[match(dist.list$Var2, env$Plot)]
  dist.list$e1 <- env[match(tolower(dist.list$Var1), tolower(rownames(env))), elev.name]
  dist.list$e2 <- env[match(tolower(dist.list$Var2), tolower(rownames(env))), elev.name]
  
  # ensure env is matching dsi.dist regarding plots/subplots
  if (is.na(dist.list$e1) || is.na(dist.list$e2)) {
    warning("The dist names not match env file plot/subplot names !\n", 
            "dist names = ", paste(firstN(unique(dist.list$Var1, dist.list$Var2), 3), collapse = ","), ", ...\n",
            "env = ", paste(firstN(rownames(env), 3), collapse = ","), ", ...\n")
    return(NA)
  }
  
  dist.list$dist <- abs(dist.list$e1-dist.list$e2)
  dist.list <- dist.list[dist.list$Var1 != dist.list$Var2, ]
  return(dist.list)
}

### Get correlation statistics ###
# One of "pearson", "kendall", or "spearman"
# analysis less than min.sample is dropped
getCorrDist <- function(dist.data, cor.method="pearson", min.sample=10, verbose=FALSE){
  corlist <- list()
  i <- 1
  for(g in unique(dist.data$gene)){
    for(m in unique(dist.data$metric)){
      if (verbose)
        cat("gene = ", g, ", metric = ", m, "\n")
      x <- dist.data[(dist.data$gene==g & dist.data$metric==m),]
      if(nrow(x) >= min.sample){
        t <- cor.test(x$value, x$dist, alternative="greater", method=cor.method)
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
