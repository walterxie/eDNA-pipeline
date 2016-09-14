# Figure 7: Correlations between spatial separation and multivariate dissimilarity index (DSI)
# Figure 8: Correlations between elevation differences and multivariate dissimilarity index (DSI)

# default by subplot
getAllDSIDist <- function(by.plot=FALSE) {
  source("R/init.R", local=TRUE)
  metrics <- c("jaccard", "horn.morisita", "bray.curtis", "beta1-1", "wt.unif", "unwt.unif")
  
  if(!exists("input.names")) stop("input names are missing !")
  output.names <- getOutputNames(input.names)
  
  all.dist.list <- list()
  for (data.id in 1:length(input.names)) {
    cm <- getIdentifiedCM(input.names[data.id], by.plot=by.plot)
    tree <- getPhyloTree(input.names[data.id])
    all.dist.list$by.plot <- by.plot
    for(metric in metrics){
      dsi.dist <- getDissimilarity(cm, tree=tree, method=metric)
      all.dist.list[[input.names[data.id]]][[metric]] <- dsi.dist
    }
  }
  return(all.dist.list) # a list of list
}


#dsi.dist.list[[input.names[data.id]]][["elevation"]] <- getElevationDistList(dsi.dist)
#dsi.dist.list[[input.names[data.id]]][["within.between"]] <- getWithinBetweenDistList(dsi.dist)
#dsi.dist.list[[input.names[data.id]]][["plot"]] <- getPlotSpatialDistList(dsi.dist)
#dsi.dist.list[[input.names[data.id]]][["subplot"]] <- getSubplotSpatialDistList(dsi.dist)

getCorSpatialDSI <- function(by.plot=TRUE, file.figure=NULL) {
  
  
  
  
}


getCorElevationDSI <- function(by.plot=TRUE, file.figure=NULL) {
  
  
}


## Get all within-between pairwise distances ###
getWithinBetweenDistDF <- function(dsi.dist){
  dist.df <- melt(as.matrix(dsi.dist))
  dist.df$p1 <- gsub("-.", "", dist.df$Var1)
  dist.df$p2 <- gsub("-.", "", dist.df$Var2)
  dist.df$s1 <- sapply(strsplit(as.character(dist.df$Var1), "-"), "[[", 2)
  dist.df$s2 <- sapply(strsplit(as.character(dist.df$Var2), "-"), "[[", 2)
  dist.df$wb <- ifelse(dist.df$p1 == dist.df$p2, "within", "between")
  dist.df <- dist.df[dist.df$Var1 != dist.df$Var2, ] # Drop same sample distances
  return(dist.df)
}

### Get elevation differences between plots/subplots ###
getElevationDistDF <- function(dsi.dist){
  dist.df <- melt(as.matrix(dsi.dist))
  #dist.df$e1 <- envdata.bysubplot$Elevation[match(dist.df$Var1, envdata.bysubplot$Plot)]
  #dist.df$e2 <- envdata.bysubplot$Elevation[match(dist.df$Var2, envdata.bysubplot$Plot)]
  dist.df$e1 <- envdata.byplot$Elevation[match(dist.df$Var1, envdata.byplot$Plot)]
  dist.df$e2 <- envdata.byplot$Elevation[match(dist.df$Var2, envdata.byplot$Plot)]
  dist.df$dist <- abs(dist.df$e1-dist.df$e2)
  dist.df <- dist.df[dist.df$Var1 != dist.df$Var2, ]
  return(dist.df)
}

### Get within-plot spatial distances ###
getSubplotSpatialDistDF <- function(dsi.dist){
  dist.df <- melt(as.matrix(dsi.dist))
  dist.df$p1 <- gsub("-.", "", dist.df$Var1)
  dist.df$p2 <- gsub("-.", "", dist.df$Var2)
  dist.df$s1 <- sapply(strsplit(as.character(dist.df$Var1), "-"), "[[", 2)
  dist.df$s2 <- sapply(strsplit(as.character(dist.df$Var2), "-"), "[[", 2)
  dist.df <- dist.df[dist.df$p1 == dist.df$p2, ] # Drop between-plot distances
  dist.df <- dist.df[dist.df$s1 != dist.df$s2, ] # Drop same sample distances
  dist.df$spair <- strSort(paste0(dist.df$s1, dist.df$s2))
  dist.df$dist <- subplot.dists$dist[match(dist.df$spair, subplot.dists$pair)]
  return(dist.df)
}

### Get between-plot spatial distances ###
getPlotSpatialDistDF <- function(dsi.dist){
  dist.df <- melt(as.matrix(dsi.dist))
  dist.df$p1 <- gsub("-.", "", dist.df$Var1)
  dist.df$p2 <- gsub("-.", "", dist.df$Var2)
  dist.df <- dist.df[dist.df$p1 != dist.df$p2, ] # Drop within-plot distances
  dist.df$pair <- plotSort(paste(dist.df$p1, dist.df$p2)) 
  dist.df$dist <- plot.dists$dist[match(dist.df$pair, plot.dists$pair)]
  return(dist.df)
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
