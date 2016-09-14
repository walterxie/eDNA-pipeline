# Figure 7: Correlations between spatial separation and multivariate dissimilarity index (DSI)
# Figure 8: Correlations between elevation differences and multivariate dissimilarity index (DSI)

# only available to subplot
getAllDSIDist <- function(save.rdata=FALSE) {
  source("R/init.R", local=TRUE)
  metrics <- c("jaccard", "horn.morisita", "bray.curtis", "beta1-1", "wt.unif", "unwt.unif")
  
  if(!exists("input.names")) stop("input names are missing !")
  output.names <- getOutputNames(input.names)
  
  by.plot=FALSE
  # calculate all dist
  cat("\nCalculate all dissimilarity ... \n")
  all.dist.list <- list()
  for (data.id in 1:length(input.names)) {
    cm <- getIdentifiedCM(input.names[data.id], by.plot=by.plot)
    tree <- getPhyloTree(input.names[data.id])
    
    for(metric in metrics){
      dsi.dist <- getDissimilarity(cm, tree=tree, method=metric)
      all.dist.list[[output.names[data.id]]][[metric]] <- dsi.dist
    }
  }
  all.dist.list$by.plot <- by.plot
  if (save.rdata)
    save(all.dist.list, file = "data/all.dist.RData")
  
  all.dist.list <- appendDistDF(output.names, metrics, all.dist.list, save.rdata=save.rdata)
  
  return(all.dist.list) # a list of list
}

appendDistDF <- function(output.names, metrics, all.dist.list=list(), save.rdata=FALSE) {
  # all within-between pairwise distances
  cat("\nCalculate all within-between pairwise distances ... \n")
  dist.data <- data.frame(check.names = F)
  for (data.id in 1:length(input.names)) {
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
  if (save.rdata)
    save(dist.data, file = "data/within.between.RData")
  
  # within-plot spatial distances
  cat("\nCalculate within-plot spatial distances ... \n")
  dist.data <- data.frame(check.names = F)
  for (data.id in 1:length(input.names)) {
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
  if (save.rdata)
    save(dist.data, file = "data/within.RData")
  
  # between-plot spatial distances
  cat("\nCalculate between-plot spatial distances ... \n")
  dist.data <- data.frame(check.names = F)
  for (data.id in 1:length(input.names)) {
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
  if (save.rdata)
    save(dist.data, file = "data/between.RData")
  
  # elevation differences between subplots
  cat("\nCalculate elevation differences between subplots ... \n")
  env.subplot <- getEnvData(by.plot=F)
  dist.data <- data.frame(check.names = F)
  for (data.id in 1:length(input.names)) {
    for(metric in metrics){
      d <- getElevationDistDF(all.dist.list[[output.names[data.id]]][[metric]], env.subplot)
      if(nrow(d) > 0){
        d$gene <- output.names[data.id]
        d$metric <- metric
      }
      dist.data <- rbind(dist.data, d)
    }
  }
  all.dist.list$elev.diff <- dist.data
  if (save.rdata)
    save(dist.data, file = "data/elev.diff.RData")
  
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


plotDistanceCorrelation <- function(dist.data, data.levels = c("16S","18S","26S","ITS","COI-300","COI-650"), 
                                    method.levels = c("Jaccard","Bray-Curtis","Morisita-Horn","beta1-1",
                                                      "unweighted Unifrac","weighted Unifrac"),
                                    x.lab="Spatial distance (m)", y.lab="Distance measure") {

  dist.data <- renderMetrics(dist.data, data.levels, method.levels)
  
  ### Distance correlation plots ###
  require(ggplot2)
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
  dat <- data.frame(x = rep(10, len), y = rep(0.5, len), vars) # x and y are coordinates relative to x/y axis scales
  dat$labs <- cor.data$cp[match(interaction(cor.data$metric, cor.data$gene), interaction(dat$metric, dat$gene))]
  
  p + geom_text(data = dat, aes(x, y, label=labs, group=NULL), size = 2)
  return(p)
}


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
  ggplot(dist.data, aes(x = wb, y = value, color = wb)) + 
    geom_boxplot(outlier.shape = 1, outlier.colour = alpha("black", 0.25)) + 
    ylab(y.lab) + xlab(x.lab) +
    facet_grid(metric ~ gene, scales = "free") +
    scale_colour_brewer(palette = "Set1") +
    theme(strip.background = element_blank(), panel.grid = element_blank(),
          legend.position = "none")
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
  
  subplot.dists  <- getSubplotDistance()
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
  
  plot.dists <- getBetweenPlotDistance()
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
