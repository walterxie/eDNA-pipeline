library(Matrix)

if (length(communityMatrix) < threshold) 
    threshold <- length(communityMatrix)

cs <- colSums(communityMatrix)

# order cm_new
ord<-order(cs, decreasing=TRUE) 
communityMatrix <- communityMatrix[,ord]    

sample_count <- length(communityMatrix[,1])
otu_count <- length(communityMatrix)
individual_count <- sum(communityMatrix)
threshold_individual_count <- sum(communityMatrix[,1:threshold])
trp <- threshold_individual_count/individual_count

sampleCounts <- apply(communityMatrix,MARGIN=2, nnzero)
