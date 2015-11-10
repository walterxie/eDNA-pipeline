library(vegan)

if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
# diss.fun = "beta1-1", "jaccard", "horn.morisita"
if(!exists("diss.fun")) diss.fun="beta1-1"
if(!exists("taxa.group")) taxa.group="all"

source("Modules/init.r")

n <- length(matrixNames)

writeDissMatrix <- function(communityMatrix, diss.fun, fname) {
  # make sure beta1-1 matrix has sample names in the same order
  communityMatrix <- communityMatrix[order(rownames(communityMatrix)),]
  
  #beta1_1 <- calculateDissimilarityMatrix(communityMatrix, printProgressBar=TRUE)
  diss.matrix <- calculateDissimilarityMatrix(communityMatrix, diss.fun)
  
  # create file for intermediate data beta1-1 matrix
  outputFile <- file.path(workingPath, "data", paste(fname, ".csv", sep = ""))
  write.cm(diss.matrix, outputFile)
} 

minOTUs <- 200
#### eDNA ####
for (expId in 1:(n-1)) {
  #### by subplot
  cat("Intermediate data: create", diss.fun, "matrix for", matrixNames[expId],
      "by subplot, rmSingleton =", rmSingleton, ", taxa.group =", taxa.group, ", otuThr =", otuThr, ".\n")
  communityMatrix <- getCommunityMatrixT(expId, FALSE, rmSingleton, taxa.group, minRow=minOTUs)
  
  fname <- paste(matrixNames[expId], postfix(taxa.group, FALSE, rmSingleton, sep="-"), diss.fun, sep = "-")
  
  if (!is.null(communityMatrix)) {
    # remove 0 row/column after merge
    communityMatrix <- prepCommunityMatrix(communityMatrix)
    
    cat(taxa.group, "subset from", matrixNames[expId], "having", ncol(communityMatrix), "OTUs", nrow(communityMatrix), "samples. \n") 
    
    writeDissMatrix(communityMatrix, diss.fun, fname)
  }
  
  #### by plot
  cat("Intermediate data: create", diss.fun, "matrix for", matrixNames[expId],
      "by plot, rmSingleton =", rmSingleton, ", taxa.group =", taxa.group, ", otuThr =", otuThr, ".\n")
  communityMatrix <- getCommunityMatrixT(expId, TRUE, rmSingleton, taxa.group, minRow=minOTUs)
  
  fname <- paste(matrixNames[expId], postfix(taxa.group, TRUE, rmSingleton, sep="-"), diss.fun, sep = "-")
  
  if (!is.null(communityMatrix)) {
    # remove 0 row/column after merge
    communityMatrix <- prepCommunityMatrix(communityMatrix)
    
    cat(taxa.group, "subset from", matrixNames[expId], "having", ncol(communityMatrix), "OTUs", nrow(communityMatrix), "samples. \n") 
    
    writeDissMatrix(communityMatrix, diss.fun, fname)
  }
} 

#### Veg ####
if (taxa.group=="all" || taxa.group=="assigned") {
  cat("Intermediate data: create", diss.fun, "matrix for", matrixNames[n], "by plot including singletons.\n")
  communityMatrix <- getCommunityMatrixT(n, TRUE, FALSE)
  
  fname <- paste(matrixNames[n], postfix("all", TRUE, FALSE, sep="-"), diss.fun, sep = "-")
  writeDissMatrix(communityMatrix, diss.fun, fname)
}
