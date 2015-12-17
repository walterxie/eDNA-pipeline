# co-occurrence network
# http://stackoverflow.com/questions/13281303/creating-co-occurrence-matrix
#
# Author: Walter Xie
# Accessed on 17 Dec 2015

#install.packages("igraph")
library(igraph)
library(reshape2)

expId=1 
isPlot=F 
rmSingleton=T
taxa.group="BACTERIA"
rankLevel="order" 
groupLevel="phylum" 

taxaAssgReads <- getTaxaAssgReads(expId, isPlot, rmSingleton, rankLevel, groupLevel, taxa.group) 

taxaAssg <- aggregate(as.formula(paste(". ~", groupLevel)), data=taxaAssgReads[,-c(1,ncol(taxaAssgReads)-1)], FUN=sum)
taxaAssg.melt <- melt(taxaAssg)

w <- dcast(taxaAssg.melt, phylum ~ variable)
x <- as.matrix(w[,-1])
x[is.na(x)] <- 0
x <- apply(x, 2,  function(x) as.numeric(x > 0))  #recode as 0/1
v <- x %*% t(x)                                   #the magic matrix 
diag(v) <- 0                                      #repalce diagonal
dimnames(v) <- list(w[, 1], w[,1])                #name the dimensions
v

g <- graph.adjacency(v, weighted=TRUE, mode ='undirected')
g <- simplify(g)
# set labels and degrees of vertices
V(g)$label <- V(g)$name
V(g)$degree <- degree(g)
g

plot(g)