
# Author: Walter Xie, Alexei Drummond
# Accessed on 9 Sep 2015

#library(vegetarian)

######## Jost diversity section #######

# communityMatrix is abundances augment in d{vegetarian}, which is transpose matrix from file
# Community data as a matrix where columns are individual species and rows are sites. 
# Matrix elements are abundance data (e.g. counts, percent cover estimates).
diversity.df <- function(communityMatrix=communityMatrix) { 
  diversity.df <- data.frame(row.names=c("gamma", "alpha", "beta"))
  
  diversity.df$'q=0' <- c(
    d(communityMatrix,lev="gamma",q=0),
    d(communityMatrix,lev="alpha",q=0),
    d(communityMatrix,lev="beta",q=0))
  
  diversity.df$'q=1' <- c(
    d(communityMatrix,lev="gamma",q=1),
    d(communityMatrix,lev="alpha",q=1),
    d(communityMatrix,lev="beta",q=1))
  
  diversity.df$'q=2' <- c(
    d(communityMatrix,lev="gamma",q=2),
    d(communityMatrix,lev="alpha",q=2),
    d(communityMatrix,lev="beta", q=2))
  
  colnames(diversity.df) <- c("$q=0$", "$q=1$", "$q=2$")
  rownames(diversity.df) <- c("$D_\\gamma(q)$", "$D_\\alpha(q)$", "$D_\\beta(q)$")
  
  return(diversity.df)
}

# abundance (reads, gamme0) per sample
# return 1-column data frame
abundancePerSample <- function(communityMatrix=communityMatrix, hasTotal) {
  if(missing(hasTotal)) hasTotal=TRUE
  
  # gamme0
  perSample <- data.frame(abundance=rowSums(communityMatrix), stringsAsFactors=FALSE)
  rownames(perSample) <- rownames(communityMatrix)
  
  if (hasTotal) {
    perSample <- rbind(perSample, sum(perSample[,1]))
    rownames(perSample)[nrow(perSample)] <- "Total"
  }
  
  return(perSample)
}

# richness (OTUs/species) per sample
# return 1-column data frame
richnessPerSample <- function(communityMatrix=communityMatrix, hasTotal) {
  if(missing(hasTotal)) hasTotal=TRUE
  
  # richness
  perSample <- data.frame(richness=rowSums(communityMatrix > 0), stringsAsFactors=FALSE)
  rownames(perSample) <- rownames(communityMatrix)
  
  if (hasTotal) {
    perSample <- rbind(perSample, sum(perSample[,1]))
    rownames(perSample)[nrow(perSample)] <- "Total"
  }
  
  return(perSample)
}

# Shannon index (gamma1) per sample
# return 1-column data frame
shannonPerSample <- function(communityMatrix=communityMatrix, digits = 2) {
  # Shannon
#  gamma1 <- function(r) d(r,lev="gamma",q=1)
#  perSample <- data.frame(Shannon=apply(communityMatrix, 1, gamma1), stringsAsFactors=FALSE)
  perSample <- data.frame(row.names=rownames(communityMatrix), stringsAsFactors=FALSE)
  for (i in 1:nrow(communityMatrix)) 
    perSample[i,1] <- d(communityMatrix[i,],lev="gamma",q=1)
  
  colnames(perSample)[1] <- "Shannon"
  
  return(round(perSample, digits))
}

######## beta 1 -1 distribution #######
# communityMatrix = t(communityMatrix), row is sample

# beta1-1 of one pair of samples r1, r2 from communityMatrix
beta1minus1.pair <- function(communityMatrix, r1, r2) {
  d(communityMatrix[c(r1, r2),],lev="beta",q=1)-1
}

# return a dist object 
beta1minus1.dist <- function(communityMatrix, printProgressBar) {    
  # including diagonal
  d.beta1minus1 <- matrix(0,nrow=nrow(communityMatrix),ncol=nrow(communityMatrix))
  colnames(d.beta1minus1) <- c(rownames(communityMatrix))
  rownames(d.beta1minus1) <- c(rownames(communityMatrix))
  # row.pairs : each row is a pair of row number of communityMatrix
  row.pairs <- t(combn(nrow(communityMatrix),2))
  
  cat("\nCalculating beta1-1 from", nrow(row.pairs), "pairs of samples.\n")
  
  if (missing(printProgressBar)) printProgressBar=nrow(row.pairs)>100
  if (printProgressBar) {
    flush.console()
    pb <- txtProgressBar(min=1, max=nrow(row.pairs), style = 3)
  }
  
  for (n in 1:nrow(row.pairs)) {
    if (printProgressBar) setTxtProgressBar(pb, n)
    d.beta1minus1[row.pairs[n,2], row.pairs[n,1]] <- d(communityMatrix[row.pairs[n,],],lev="beta",q=1)-1
  }
  if (printProgressBar) close(pb)
  
  return (as.dist(d.beta1minus1))
}


