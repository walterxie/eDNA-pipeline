# OPTIMAL NUMBER OF CLUSTERS  
# http://enterotype.embl.de/enterotypes_tutorial.sanger.R
#
# Author: Walter Xie
# Accessed on 23 Nov 2015

#install.packages("cluster")
#install.packages("clusterSim")
library(cluster)
library(clusterSim)
library(ade4)
library(ggplot2)
library(ellipse)

# Jensen-Shannon divergence (JSD) (Endres & Schindelin, 2003)
# Kullback-Leibler divergence
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}
#data.dist=dist.JSD(data)

# Partitioning around medoids (PAM) clustering algorithm to cluster the abundance profiles. 
# PAM derives from the basic k-means algorithm, but has the advantage that it supports any 
# arbitrary distance measure and is more robust than k-means. It is a supervised procedure, 
# where the predetermined number of clusters is given as input to the procedure, which then 
# partitions the data into exactly that many clusters.
pam.clustering=function(x,k, as.vector=TRUE) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster=pam(as.dist(x), k, diss=TRUE)
  if (as.vector) {
    return(as.vector(cluster$clustering))
  }
  return(cluster)
}
#data.cluster=pam.clustering(data.dist, k=3)

# n.max: max number of clusters to run
run.pam.cl=function(data, data.dist, n.max=20) { 
  require(clusterSim)
  #nclusters = index.G1(t(data), data.cluster, d = data.dist, centrotypes = "medoids")
  nclusters=NULL
  for (k in 1:n.max) { 
    if (k==1) {
      nclusters[k]=NA 
    } else {
      data.cluster_temp=pam.clustering(data.dist, k)
      nclusters[k]=index.G1(t(data), data.cluster_temp, d=data.dist, centrotypes = "medoids")
    }
  }
  return(nclusters)
}
#nclusters=run.pam.cl(data, data.dist, n.max)

plotOptClusters <- function(figPath, nclusters, n.max=20) {
  # the optimal number of clusters
  pdf(file.path(figPath, "optimal-clusters-JSD.pdf"), width=6, height=6)
  plot(nclusters, type="h", xlab="number of clusters", ylab="CH index",
       main=paste0("Optimal number of clusters (Jensen-Shannon divergence)"), xaxt="n")
  axis(1, at = seq(1, n.max, by = 2))
  invisible(dev.off())   
}

plotClusterAbundence <- function(figPath, data, data.cluster) {
  k=max(data.cluster)
  da.cl <- as.data.frame(cbind(t(data), cluster=data.cluster))
  da.cl.melt <- melt(da.cl, id=c("cluster"))
  da.cl.melt <- da.cl.melt[da.cl.melt$value>0,]
  
  for (cl in sort(unique(data.cluster))) {
    gg <- da.cl.melt[da.cl.melt$cluster==cl,]
    
    cat("Cluster", cl, "smaples:", paste(rownames(da.cl[da.cl$cluster==cl,]), collapse = ", "), ".\n") 
    
    p <- ggplot(gg, aes(x= reorder(variable, value, median, order=TRUE), y=value)) + 
      geom_boxplot(aes(fill = factor(variable))) + 
      ggtitle(paste("Cluster", cl)) + xlab("") + ylab("Abundence") + 
      theme_bw() + guides(fill=FALSE) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) 
    
    pdf(file.path(figPath, paste0("enterotype-",k,"-",cl,".pdf")), width=8, height=6)
    print(p)
    invisible(dev.off()) 
  }
}

plotEnterotypes <- function(figPath, data, data.dist, data.cluster, percent=0.01, nf=3, addLabel=TRUE) {
  k=max(data.cluster)
  
  if (percent > 0)
    data=noise.removal(data, percent)
  
  plotClusterAbundence(figPath, data, data.cluster)
  
  #Between-class analysis (BCA) was performed to support the clustering and identify the drivers for the enterotypes. 
  obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=nf)
  obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=nf) 
  cat("\nplot", k, "clusters between-class analysis (BCA).\n") 
  
  plotScatterEllipse(as.data.frame(obs.bet$ls), data.cluster, "Between class analysis", 
                     file.path(figPath, paste0("enterotypes-between-class-",k,".pdf")), addLabel)

  #principal coordinates analysis (PCoA) of a Euclidean distance matrix
  obs.pcoa=dudi.pco(data.dist, scannf=F, nf=nf)
  cat("\nplot", k, "clusters principal coordinates analysis (PCoA).\n") 
  
  plotScatterEllipse(obs.pcoa$li, data.cluster, "Principal coordinates analysis", 
                     file.path(figPath, paste0("enterotypes-principal-coordiante-",k,".pdf")), addLabel)
}

# http://stackoverflow.com/questions/2397097/how-can-a-data-ellipse-be-superimposed-on-a-ggplot2-scatterplot
plotScatterEllipse <- function(data.obs, data.cluster, title, figPathName, addLabel=TRUE) {
  gg <- cbind(data.obs, cluster=data.cluster)
  colnames(gg)[1:2] <- c("PC1", "PC2")
#  gg$species <- paste(sapply(strsplit(rownames(gg), "_"), "[[", 1), sapply(strsplit(rownames(gg), "_"), "[[", 2), sep=".")
  gg$cluster <- factor(gg$cluster, levels = sort(unique(gg$cluster)))
  
  p <- ggplot(gg, aes(x=PC1, y=PC2, color=factor(cluster))) + 
    geom_point(size=3) + #aes(shape=factor(species)),
    scale_colour_brewer(name="cluster", palette="Set1") + #scale_fill_manual(name = "cluster", values = myPalette)  +
#    scale_shape_manual(name ="species", values=1:length(unique(gg$species))) +
    stat_ellipse(type = "t", linetype = 2) +
    geom_hline(yintercept=0,linetype=2) + 
    geom_vline(xintercept=0,linetype=2) +
    ggtitle(title) +
    theme_bw() + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), panel.background = element_blank()) 

  if (addLabel) {
    # Error in row.names(gg) : object 'gg' not found
    p <- p + geom_text(aes(color=factor(cluster), label=rownames(gg)),size = 3, hjust=-0.1, vjust = -0.2, alpha = 0.5) 
  }
  
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  
  pdf(figPathName, width=8, height=6)
  print(grid.draw(gt))
  invisible(dev.off()) 
}

validateOptClusters <- function(data, data.dist, data.cluster) {
  k=max(data.cluster)
   
  nclusters = index.G1(t(data), data.cluster, d = data.dist, centrotypes = "medoids")
  
  # Observations with a large s(i) (almost 1) are very well clustered, a small s(i) (around 0) means 
  # that the observation lies between two clusters, and observations with a negative s(i) are probably 
  # placed in the wrong cluster.
  obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
  cat("obs.silhouette =", obs.silhouette, ", when number of clusters =", k, "\n") 
}

# advise to apply this function to data generated using short sequencing technologes, like Illumina or Solid
noise.removal <- function(Matrix, percent=0.01){
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  cat("percent = ", percent, ", remove", nrow(Matrix)-nrow(Matrix_1), "rows.\n")
  return(Matrix_1)
}

######### vegan cascadeKM ##########
plotCascadeKM <- function(figPath, data, n.max) {
  library(vegan)
  
  ccas <- cascadeKM(t(data), 2, n.max)
  #ccas
  pdf(file.path(figPath, paste0("optimal-clusters-cascadeKM.pdf")), width=6, height=6)
  plot(ccas, sortq=TRUE)
  invisible(dev.off()) 
  
  ccas <- cascadeKM(t(data), 2, 20, criterion = "ssi")
  pdf(file.path(figPath, paste0("optimal-clusters-cascadeKM-ssi.pdf")), width=6, height=6)
  plot(ccas, sortq=TRUE)
  invisible(dev.off())
}

######### deprecated ##########
plotEnterotypes.bak <- function(figPath, data, data.dist, data.cluster, k, percent=0.01) {
  if (percent > 0)
    data=noise.removal(data, percent)
  
  ## plot 1
  obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
  obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k) 
  
  pdf(file.path(figPath, "enterotypes-between-class.pdf"), width=6, height=6)
  #  dev.new()
  s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F,sub="Between-class analysis")
  invisible(dev.off()) 
  
  #plot 2
  obs.pcoa=dudi.pco(data.dist, scannf=F, nf=k)
  
  pdf(file.path(figPath, "enterotypes-principal-coordiante.pdf"), width=6, height=6)
  #  dev.new()
  s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F,sub="Principal coordiante analysis")
  invisible(dev.off())
}
