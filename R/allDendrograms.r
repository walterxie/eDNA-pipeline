
if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("diss.fun")) diss.fun="beta1-1"
if(!exists("taxa.group")) taxa.group="all"

n <- length(matrixNames) 

source("Modules/init.R", local=TRUE)

# function to get color labels
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- mypalette[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}

n_cluster=9
mypalette <- rainbow(n_cluster)

cat("Analysis:", n_cluster, "-cluster UPGMA dendrogram of", taxa.group, "taxa group(s) using", diss.fun, "dissimilarity. \n")

######## dendrograms #######	
for (expId in 1:n) {	
  # isPlot, rmSingleton, taxa.group, are fixed in init, when expId == n
  diss <- getDissimilarityMatrix(expId, TRUE, rmSingleton, diss.fun, taxa.group)

  if (expId == n) {
    fname <- paste("dendrogram", matrixNames[expId], postfix("all", TRUE, FALSE, sep="-"), diss.fun, sep = "-")
  } else {
    fname <- paste("dendrogram", matrixNames[expId], postfix(taxa.group, TRUE, rmSingleton, sep="-"), diss.fun, sep = "-")
  }
  
  hc <- hclust(as.dist(diss), "ave") # "average" (= UPGMA)
  
  # cut dendrogram in n_cluster clusters
  clusMember = cutree(hc, n_cluster)
  # using dendrogram objects
  hcd = as.dendrogram(hc)
  # using dendrapply
  clusDendro = dendrapply(hcd, colLab)
  
  pdf(paste(workingPath, figDir, "/", fname, ".pdf", sep = ""), width=6, height=6)
  
  par(mar=c(5,3,3,1)) 
  plot(clusDendro, xlab="", sub ="", main=paste(n_cluster, "-cluster UPGMA dendrogram of", matrixNames[expId], sep=" "))
  
  invisible(dev.off())
}
