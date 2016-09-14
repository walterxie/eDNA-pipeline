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
  
  return(all.dist.list) # a list of list
}


getCorSpatialDSI <- function(by.plot=TRUE, file.figure=NULL) {
  
  
  
  
}


getCorElevationDSI <- function(by.plot=TRUE, file.figure=NULL) {
  
  
}


## Get all within-between pairwise distances ###
getWithinBetweenDistDF <- function(dsi.dist){
  dist.list <- melt(as.matrix(dsi.dist))
  dist.list$p1 <- gsub("-.", "", dist.list$Var1)
  dist.list$p2 <- gsub("-.", "", dist.list$Var2)
  dist.list$s1 <- sapply(strsplit(as.character(dist.list$Var1), "-"), "[[", 2)
  dist.list$s2 <- sapply(strsplit(as.character(dist.list$Var2), "-"), "[[", 2)
  dist.list$wb <- ifelse(dist.list$p1 == dist.list$p2, "within", "between")
  dist.list <- dist.list[dist.list$Var1 != dist.list$Var2, ] # Drop same sample distances
  return(dist.list)
}

### Get elevation differences between plots/subplots ###
getElevationDistDF <- function(dsi.dist, env, env.names=c("Plot","Elevation")){
  if (! any(env.names %in% colnames(env)) )
    stop("Cannot find column names from env file ! ", paste(env.names, collapse = ","))
  dist.list <- melt(as.matrix(dsi.dist))
  #dist.list$e1 <- env$Elevation[match(dist.list$Var1, env$Plot)]
  #dist.list$e2 <- env$Elevation[match(dist.list$Var2, env$Plot)]
  dist.list$e1 <- env[match(dist.list$Var1, env[, env.names[1]]), env.names[2]]
  dist.list$e2 <- env[match(dist.list$Var2, env[, env.names[1]]), env.names[2]]
  dist.list$dist <- abs(dist.list$e1-dist.list$e2)
  dist.list <- dist.list[dist.list$Var1 != dist.list$Var2, ]
  return(dist.list)
}

### Get within-plot spatial distances ###
getWithinDistDF <- function(dsi.dist){
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
  dist.list <- melt(as.matrix(dsi.dist))
  dist.list$p1 <- gsub("-.", "", dist.list$Var1)
  dist.list$p2 <- gsub("-.", "", dist.list$Var2)
  dist.list <- dist.list[dist.list$p1 != dist.list$p2, ] # Drop within-plot distances
  dist.list$pair <- plotSort(paste(dist.list$p1, dist.list$p2)) 
  
  plot.dists <- getBetweenPlotDistance()
  dist.list$dist <- plot.dists$dist[match(dist.list$pair, plot.dists$pair)]
  return(dist.list)
}

### Get correlation statistics ###
get.corr <- function(dist.data){
  corlist <- list()
  i <- 1
  for(g in unique(dist.data$gene)){
    for(m in unique(dist.data$metric)){
      print(paste(g, m))
      x <- dist.data[(dist.data$gene==g & dist.data$metric==m),]
      if(nrow(x) > 10){
        t <- cor.test(x$value, x$dist, alternative="greater", method=c("pearson"))
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
