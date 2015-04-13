## rotate cm for library(untb), and remove all colSums(communityMatrix) == 0, length >= threshold

communityMatrix <- as.data.frame(t(communityMatrix)) # rotate  

names(communityMatrix) <- sub(" ", ".", names(communityMatrix))  
if (length(communityMatrix) < threshold) 
    threshold <- length(communityMatrix)

communityMatrix <- communityMatrix[,colSums(communityMatrix) > 0] # remove all taxa have no values (all 0) 

print(paste("transposed community matrix ", matrixName, ": columns = ", ncol(communityMatrix), ", rows = ", nrow(communityMatrix), sep=""))   

