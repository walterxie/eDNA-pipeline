
# Author: Walter Xie, Alexei Drummond
# Accessed on 9 Sep 2015

library(vegetarian)

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
abundancePerSample <- function(communityMatrix, hasTotal) {
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
richnessPerSample <- function(communityMatrix, hasTotal) {
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
shannonPerSample <- function(communityMatrix, digits = 2) {
  # Shannon
  #  gamma1 <- function(r) d(r,lev="gamma",q=1)
  #  perSample <- data.frame(Shannon=apply(communityMatrix, 1, gamma1), stringsAsFactors=FALSE)
  perSample <- data.frame(row.names=rownames(communityMatrix), stringsAsFactors=FALSE)
  for (i in 1:nrow(communityMatrix)) 
    perSample[i,1] <- d(communityMatrix[i,],lev="gamma",q=1)
  
  colnames(perSample)[1] <- "Shannon"
  
  return(round(perSample, digits))
}

######## alpha1 #######
alpha1 <- function(communityMatrix) {    
  # including diagonal
  m.alpha1 <- matrix(0,nrow=nrow(communityMatrix),ncol=1)	
  rownames(m.alpha1) <- c(rownames(communityMatrix))
  for(i in 1:nrow(communityMatrix)){				
    m.alpha1[i,1] <- d(communityMatrix[i,],lev="gamma",q=1)				
  }
  
  return (m.alpha1) # 1 col matrix
}

######## Pair-wise turnovers #######
# communityMatrix = t(communityMatrix), row is sample

# return a matrix cols and rows are sample names 
# printProgressBar = TRUE/FALSE, if missing, =nrow(row.pairs)>100
calculateDissimilarityMatrix <- function(communityMatrix, diss.fun="beta1-1", printProgressBar) {    
  # including diagonal
  diss.matrix <- matrix(0,nrow=nrow(communityMatrix),ncol=nrow(communityMatrix))
  colnames(diss.matrix) <- c(rownames(communityMatrix))
  rownames(diss.matrix) <- c(rownames(communityMatrix))
  # row.pairs : each row is a pair of row number of communityMatrix
  row.pairs <- t(combn(nrow(communityMatrix),2))
  
  cat("\nCalculating", diss.fun, "from", nrow(row.pairs), "pairs of samples.\n")
  
  if (missing(printProgressBar)) printProgressBar=nrow(row.pairs)>100
  if (printProgressBar) {
    flush.console()
    pb <- txtProgressBar(min=1, max=nrow(row.pairs), style = 3)
  }
  
  for (n in 1:nrow(row.pairs)) {
    if (printProgressBar) setTxtProgressBar(pb, n)
    if (diss.fun=="jaccard") {
      # Jaccard
      diss.matrix[row.pairs[n,2], row.pairs[n,1]] <- vegdist(communityMatrix[row.pairs[n,],], method="jaccard", binary=TRUE)
    } else if (diss.fun=="horn.morisita") {
      # Horn-Morisita
      diss.matrix[row.pairs[n,2], row.pairs[n,1]] <- vegdist(communityMatrix[row.pairs[n,],], method="horn", binary=FALSE)
    } else if (diss.fun=="bray.curtis") {
      # Bray-Curtis
      diss.matrix[row.pairs[n,2], row.pairs[n,1]] <- vegdist(communityMatrix[row.pairs[n,],])
    } else { # diss.fun="beta1-1"
      # beta1-1
      diss.matrix[row.pairs[n,2], row.pairs[n,1]] <- d(communityMatrix[row.pairs[n,],],lev="beta",q=1)-1
    }
  }
  if (printProgressBar) close(pb)
  
  return (diss.matrix)
}

# Returns a distance matrix composed of pair-wise turnovers
# Input: a community matrix
# Depends on: vegetarian library 
TurnoverDist<-function(input.table){   
  to.table<-matrix(0,nrow=nrow(input.table),ncol=nrow(input.table))
  for(i in 1:nrow(input.table)){
    for(j in 1:nrow(input.table)){
      # For numerous communities of equal weights, the numbers equivalent of 
      # the Shannon beta diversity and the number of samples (N) can be used to 
      # calculate the turnover rate per sample (Equation 25 from Jost 2007, Harrison et al. 1992)
      to.table[i,j]<-turnover(input.table[c(i,j),])
    }
  }
  
  d <- as.dist(to.table)
  attr(d, "Labels") <- dimnames(input.table)[[1L]]
  
  return(d)
}


# COMPUTE TURNOVER TABLE
#
#d.turnover <- TurnoverDist(communityMatrix)
#

# COMPUTER HORN-MORISITA OVERLAPS

#library(vegan)
#d.hornMorisita <- vegdist(communityMatrix, method="horn", binary=FALSE)
#d.brayBin <- vegdist(communityMatrix, method="bray", binary=TRUE)


#library(untb)
#cm_counts <- count(colSums(communityMatrix))
#theta <- round(optimal.theta(cm_counts),2)
