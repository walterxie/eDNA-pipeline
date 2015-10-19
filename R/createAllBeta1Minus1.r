library(vegan)

if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97

source("Modules/init.r")

cat("Intermediate data: create all beta1-1 matrix, rmSingleton =", rmSingleton, ", otuThr =", otuThr, "\n") 

n <- length(matrixNames)
# main
for (expId in 1:n) {
  min2 <- rmSingleton
  if (expId == n) 
    min2 <- FALSE
  
  #### by subplot
  if (expId < n) {
    cat("Create beta1-1 matrix by subplot for", matrixNames[expId],
        ", rmSingleton =", min2, ", otuThr =", otuThr, ".\n")
    communityMatrix <- getCommunityMatrixT(expId, FALSE, min2)
    
    #beta1_1 <- beta1minus1(communityMatrix, printProgressBar=TRUE)
    beta1_1 <- beta1minus1(communityMatrix)
    
    # create file for intermediate data beta1-1 matrix
    outputFile <- paste(workingPath, "data/", postfix(matrixNames[expId], FALSE, min2, sep="-"), "-beta1-1.csv", sep = "")
    write.cm(beta1_1, outputFile)
  }
  
  #### by plot
  cat("Create beta1-1 matrix by plot for", matrixNames[expId],
      ", rmSingleton =", min2, ", otuThr =", otuThr, ".\n")
  communityMatrix <- getCommunityMatrixT(expId, TRUE, min2)
  
  #beta1_1 <- beta1minus1(communityMatrix, printProgressBar=TRUE)
  beta1_1 <- beta1minus1(communityMatrix)
  
  # create file for intermediate data beta1-1 matrix
  outputFile <- paste(workingPath, "data/", postfix(matrixNames[expId], TRUE, min2, sep="-"), "-beta1-1.csv", sep = "")
  write.cm(beta1_1, outputFile)
} 



