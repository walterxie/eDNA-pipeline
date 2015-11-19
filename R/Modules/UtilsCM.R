# Author: Walter Xie, Andrew Dopheide
# Accessed on 19 Nov 2015


# remove singletons from rows/cols
rmSingletons <- function(communityMatrix, MARGIN=1, verbose=TRUE) {
  if (is.element(1, MARGIN)) {
    singletons <- which(rowSums(communityMatrix)==1)
    communityMatrix <- communityMatrix[-singletons,]
    msg <- "from rows"
  } else if (is.element(2, MARGIN)) {
    singletons <- which(colSums(communityMatrix)==1)
    communityMatrix <- communityMatrix[,-singletons]
    msg <- "from columns"
  }
  
  if(verbose) 
    cat("Remove", length(singletons) ,"singletons", msg, "!\n")
  
  communityMatrix
}

# communityMatrix: data frame
transposeCM <- function(communityMatrix=communityMatrix) {
  if (!all(sapply(communityMatrix, is.numeric))) 
    stop("All community matrix elements have to be numeric type") 
  
  communityMatrixT <- as.data.frame(t(as.matrix(communityMatrix))) # transpose  
}



### Combine columns by sample name
mergeCMSample <- function(communityMatrix=communityMatrix, sep) {
  if(missing(sep)) sep="-"
  
  colnames(communityMatrix) <- sapply(strsplit(colnames(communityMatrix), sep), "[[", 1) # Strip subplot letter by sep
  communityMatrix1 <- data.frame(matrix(ncol = 0, nrow = nrow(communityMatrix))) # Empty data.frame with required number of rows
  for(col in unique(colnames(communityMatrix))){
    cols <- communityMatrix[grep(col, colnames(communityMatrix))] # Find each pair of subplot columns
    cols1 <- as.data.frame(rowSums(cols)) # Add subplot pair together
    colnames(cols1) <- col
    communityMatrix1 <- cbind(communityMatrix1, cols1)
  }
  
  return(communityMatrix1)
}

# rowThr, colThr: remove row or/and column sum <= rowThr, colThr in communityMatrix,
# if rowThr, colThr = 0, then remove empty rows or/and columns
# mostAbundThr: keep the most abundant (mostAbundThr) OTUs only 
# return preprocessed communityMatrix
preprocessCM <- function(communityMatrix, keepSingleton, rowThr, colThr, mostAbundThr) { 
  # keepSingleton, abundPercThr, transverse, verbose
  # args <- list(...)
  if(missing(keepSingleton)) keepSingleton=TRUE
  if(missing(rowThr)) rowThr=0
  if(missing(colThr)) colThr=0
  if(missing(mostAbundThr)) mostAbundThr=0
  
  cat("keepSingleton = ", keepSingleton, "; rowThr = ", rowThr, "; colThr = ", colThr, 
      "; mostAbundThr = ", mostAbundThr, "\n") 
  cat("Original community matrix : samples = ", ncol(communityMatrix), ", OTUs/taxa = ", nrow(communityMatrix), ".\n") 
  
  # singletons
  if (!keepSingleton) {
    singletons <- which(rowSums(communityMatrix)==1)
    communityMatrix <- communityMatrix[-singletons,]
    cat("Remove", length(singletons) ,"singletons !\n")
    rm(singletons)		
  }	
  
  # this must be in front of filter column/row to avoid empty column
  if (mostAbundThr > 0) {
    communityMatrix <- keepMostAbundantRows(communityMatrix, mostAbundThr=mostAbundThr)
  }  
  
  communityMatrix <- rmVectorFromCM(communityMatrix, vectorThr=colThr, MARGIN=2)
  # filter column first to avoid empty rows after columns remvoed
  communityMatrix <- rmVectorFromCM(communityMatrix, vectorThr=rowThr, MARGIN=1)
  
  # summary
  cat("Processed community matrix : samples = ", ncol(communityMatrix), ", OTUs/taxa = ", nrow(communityMatrix), ".\n") 
  
  return(communityMatrix)
}


# MARGIN: 1 indicates rows, 2 indicates columns, c(1, 2) indicates rows and columns.
# remove row or/and column sum <= vectorThr in communityMatrix
rmVectorFromCM <- function(communityMatrix, vectorThr, MARGIN) {
  if(missing(vectorThr)) vectorThr=0
  if(missing(MARGIN)) MARGIN=c(1,2)
  
  row.count <- nrow(communityMatrix)	
  column.count <- ncol(communityMatrix)
  
  # 1 row/col issue
  if (row.count > 2 && column.count > 2) {
    if (is.element(1, MARGIN)) {
      communityMatrix <- communityMatrix[rowSums(communityMatrix) > vectorThr,]
    } else if (is.element(2, MARGIN)) {
      communityMatrix <- communityMatrix[,colSums(communityMatrix) > vectorThr]
    }
    
    if (row.count != nrow(communityMatrix)) 
      cat("Warning : remove", row.count-nrow(communityMatrix), "row(s) whose row sum <=", vectorThr, "!\n")
    if (column.count != ncol(communityMatrix)) 
      cat("Warning : remove", column.count-ncol(communityMatrix), "column(s) whose column sum <=", vectorThr, "!\n")
  } else {
    cat("Warning : skip filtering for 1 row/column data frame !\n")
  }
  
  return(communityMatrix)
}


