# OPTIMAL NUMBER OF CLUSTERS  
# http://enterotype.embl.de/enterotypes_tutorial.sanger.R
#
# Author: Walter Xie
# Accessed on 23 Nov 2015

#install.packages("cluster")
#install.packages("clusterSim")
library(cluster)
library(clusterSim)

# change path here
workingPath <- "~/WorkSpace/eDNA-pipeline/Vegetation"
setwd(workingPath)


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

# Partitioning around medoids (PAM) clustering algorithm to cluster the abundance profiles. 
# PAM derives from the basic k-means algorithm, but has the advantage that it supports any 
# arbitrary distance measure and is more robust than k-means. It is a supervised procedure, 
# where the predetermined number of clusters is given as input to the procedure, which then 
# partitions the data into exactly that many clusters.
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
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

upgma.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

run.upgma.cl=function(data, data.dist, n.max=20) { 
  require(clusterSim)
  #nclusters = index.G1(t(data), data.cluster, d = data.dist, centrotypes = "medoids")
  nclusters=NULL
  for (k in 1:n.max) { 
    if (k==1) {
      nclusters[k]=NA 
    } else {
      data.cluster_temp=hclust(data.dist, "ave") # "average" (= UPGMA)
      # cut dendrogram in n_cluster clusters
      clusMember = cutree(data.cluster_temp, k)
      nclusters[k]=index.G1(t(data), data.cluster_temp, d=data.dist, centrotypes = "medoids")
    }
  }
  return(nclusters)
}

source("../R/Modules/Diversities.R")

###### load SBA 
sbaFile <- file.path(workingPath, "LBI_Trees_Saplings_SBA.csv")
sba <- read.csv(sbaFile, head=TRUE, check.names=FALSE, row.names = 1, stringsAsFactors=FALSE)
cat("\nUpload SBA file: ", nrow(sba), "rows", ncol(sba), "columns from", sbaFile, "\n")
#sba

n.max=20
###### Jensen-Shannon Divergence
dist.fun="Jensen-Shannon divergence"
sba.dist=dist.JSD(sba)

nclusters=run.pam.cl(sba, sba.dist, n.max)

# the optimal number of clusters
pdf(file.path(workingPath, "optimal-clusters-JSD.pdf"), width=6, height=6)
plot(nclusters, type="h", xlab="number of clusters", ylab="CH index",
     main=paste0("Optimal number of clusters (", dist.fun, ")"), xaxt="n")
axis(1, at = seq(1, n.max, by = 2))
invisible(dev.off())   

k=2
sba.cluster=pam.clustering(sba.dist, k=k)
# Observations with a large s(i) (almost 1) are very well clustered, a small s(i) (around 0) means 
# that the observation lies between two clusters, and observations with a negative s(i) are probably 
# placed in the wrong cluster.
obs.silhouette=mean(silhouette(sba.cluster, sba.dist)[,3])
cat("obs.silhouette =", obs.silhouette, ", when optimal number of clusters =", k, "\n") 

###### beta1-1
dist.fun="beta1-1"
diss.matrix <- calculateDissimilarityMatrix(t(sba), diss.fun=dist.fun)
diss.matrix <- diss.matrix * 10
sba.dist=as.dist(diss.matrix)

nclusters=run.upgma.cl(sba, sba.dist, n.max)

# the optimal number of clusters
pdf(file.path(workingPath, paste0("optimal-clusters-", dist.fun, ".pdf")), width=6, height=6)
plot(nclusters, type="h", xlab="number of clusters", ylab="CH index",
     main=paste0("Optimal number of clusters (", dist.fun, ")"), xaxt="n")
axis(1, at = seq(1, n.max, by = 2))
invisible(dev.off())   


###### jaccard
dist.fun="jaccard"
diss.matrix <- calculateDissimilarityMatrix(t(sba), diss.fun=dist.fun)
sba.dist=as.dist(diss.matrix)

nclusters=run.upgma.cl(sba, sba.dist, n.max)

# the optimal number of clusters
pdf(file.path(workingPath, paste0("optimal-clusters-", dist.fun, ".pdf")), width=6, height=6)
plot(nclusters, type="h", xlab="number of clusters", ylab="CH index",
     main=paste0("Optimal number of clusters (", dist.fun, ")"), xaxt="n")
axis(1, at = seq(1, n.max, by = 2))
invisible(dev.off())   

######### vegan cascadeKM ##########
library(vegan)

ccas <- cascadeKM(t(sba), 2, n.max)
#ccas
pdf(file.path(workingPath, paste0("optimal-clusters-cascadeKM.pdf")), width=6, height=6)
plot(ccas, sortq=TRUE)
invisible(dev.off()) 

