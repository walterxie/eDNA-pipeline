library(vegan)

if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97

if(!exists("diss.fun")) diss.fun="beta1-1"
if(!exists("taxa.group")) taxa.group="all"

source("Modules/init.r")

n <- length(matrixNames)

#### eDNA ####
for (expId in 1:(n-1)) {
  #### by subplot
  cat("Intermediate data: create", diss.fun, "matrix for", matrixNames[expId],
      "by subplot, rmSingleton =", rmSingleton, ", taxa.group =", taxa.group, ", otuThr =", otuThr, ".\n")
  communityMatrix <- getCommunityMatrixT(expId, FALSE, rmSingleton, taxa.group)
  
  # make sure beta1-1 matrix has sample names in the same order
  communityMatrix <- communityMatrix[order(rownames(communityMatrix)),]
  
  #beta1_1 <- calculateDissimilarityMatrix(communityMatrix, printProgressBar=TRUE)
  diss.matrix <- calculateDissimilarityMatrix(communityMatrix, diss.fun)
  
  fname <- paste(matrixNames[expId], postfix(taxa.group, FALSE, rmSingleton, sep="-"), diss.fun, sep = "-")
  # create file for intermediate data beta1-1 matrix
  outputFile <- paste(workingPath, "data/", fname, ".csv", sep = "")
  write.cm(diss.matrix, outputFile)
  
  #### by plot
  cat("Intermediate data: create", diss.fun, "matrix for", matrixNames[expId],
      "by plot, rmSingleton =", rmSingleton, ", taxa.group =", taxa.group, ", otuThr =", otuThr, ".\n")
  communityMatrix <- getCommunityMatrixT(expId, TRUE, rmSingleton, taxa.group)
  
  #beta1_1 <- calculateDissimilarityMatrix(communityMatrix, printProgressBar=TRUE)
  diss.matrix <- calculateDissimilarityMatrix(communityMatrix, diss.fun)
  
  fname <- paste(matrixNames[expId], postfix(taxa.group, TRUE, rmSingleton, sep="-"), diss.fun, sep = "-")
  # create file for intermediate data beta1-1 matrix
  outputFile <- paste(workingPath, "data/", fname, ".csv", sep = "")
  write.cm(diss.matrix, outputFile)
} 

#### Veg ####
cat("Intermediate data: create", diss.fun, "matrix for", matrixNames[n], "by plot including singletons.\n")
communityMatrix <- getCommunityMatrixT(expId, TRUE, FALSE)

#beta1_1 <- calculateDissimilarityMatrix(communityMatrix, printProgressBar=TRUE)
diss.matrix <- calculateDissimilarityMatrix(communityMatrix, diss.fun)

fname <- paste(postfix(matrixNames[n], TRUE, FALSE, sep="-"), diss.fun, sep = "-")
# create file for intermediate data beta1-1 matrix
outputFile <- paste(workingPath, "data/", fname, ".csv", sep = "")
write.cm(diss.matrix, outputFile)

