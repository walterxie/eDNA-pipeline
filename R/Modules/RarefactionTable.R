
library(vegan)


if(!exists("matrixNames")) stop("matrix names are missing !")

# return 
#  alpha0 alpha1
#CM30c28     ?     ?
rrarefyPerSample <- function(communityMatrix, sampleSize, reps) {
  n <- nrow(communityMatrix)
  rdMeanTable <- data.frame(row.names=rownames(communityMatrix), alpha0=rep(NA,n), alpha1=rep(NA,n), check.names=FALSE)
  for (r in 1:n) {
    if (sum(communityMatrix[r,]) >= sampleSize) {
      rdt0 <- c()
      rdt1 <- c()
      for (i in 1:reps) {
        cmrep <- rrarefy(communityMatrix[r,], sampleSize)
        rd <- d(cmrep,lev="gamma",q=0)
        rdt0 <- c(rdt0, rd)
        rd <- d(cmrep,lev="gamma",q=1)
        rdt1 <- c(rdt1, rd)
      }
      rdMeanTable[r,"alpha0"] <- mean(rdt0)
      #rdSterr <- sd(rdt)/sqrt(reps)
      rdMeanTable[r,"alpha1"] <- mean(rdt1)
    }
  }
  rdMeanTable
}

createRDPerSampleTable <- function(expId, isPlot, rmSingleton, taxa.group, pathFileStem, reps=10, min.size=100, verbose=T) {
  communityMatrix <- getCommunityMatrixT(expId, isPlot, rmSingleton, taxa.group)
  
  rare.max <- max(rowSums(communityMatrix))
  rare.min <- min(rowSums(communityMatrix))
  
  if (rare.min < min.size) {
    cat("\nWarning: min sample size", rare.min, "<", min.size, ", skip", taxa.group, "subset from", matrixNames[expId], ".\n")
    return(FALSE)
  }
  sample.sizes <- c(round(exp(seq(log(1), log(rare.min), length.out = 6)), digits = 0)[1:5], 
                    round(exp(seq(log(rare.min), log(rare.max), length.out = 6)), digits = 0))
  
  if (verbose) {
    cat("Subsampling: reps =", reps, ", min size allowed", min.size,".\n") 
    cat("Rarefaction :", matrixNames[expId], taxa.group, ", min =", rare.min, ", max =", rare.max, ".\n") 
  }
  
  alpha0MeanTable <- data.frame(row.names=rownames(communityMatrix), check.names=FALSE)
  alpha1MeanTable <- data.frame(row.names=rownames(communityMatrix), check.names=FALSE)
  
  for (ss in sample.sizes) {
    if (verbose)
      cat("sample size =", ss, ".\n") 
    
    rdMeanTable <- rrarefyPerSample(communityMatrix, ss, reps)
    alpha0MeanTable[,paste("size.", ss, sep="")] <- rdMeanTable$alpha0
    alpha1MeanTable[,paste("size.", ss, sep="")] <- rdMeanTable$alpha1
  }
  
  outputFile <- paste(pathFileStem, "rare", "alpha0", "table.csv", sep = "-") 
  write.csv(alpha0MeanTable, outputFile, row.names=TRUE, quote=FALSE)
  
  if (verbose) 
    cat("Write file", outputFile, ".\n")

  outputFile <- paste(pathFileStem, "rare", "alpha1", "table.csv", sep = "-")
  write.csv(alpha1MeanTable, outputFile, row.names=TRUE, quote=FALSE)
  
  if (verbose) 
    cat("Write file", outputFile, ".\n")
  
  return(TRUE)
}

######## old code #######
levels = rep(c("gamma","alpha","beta"),3)
qs = rep(0:2,each=3)

getMinSizeAllSites <- function(communityMatrix) {
  minSample <- min(99999, min(rowSums(communityMatrix)))  		
  
  sampleSize <- minSample#floor(minSample/100) * 100 # make sure subsampling run  
  cat("Community matrix: min sample size per site =", minSample, ", take sample size per site =", sampleSize, "\n") 
  
  return (sampleSize)
}

# 9 rows of Jost diversities, columns are mean, sterrD   
rdiversityTable <- function(communityMatrix, sampleSize, reps) {
  size = length(levels)
  
  rdt <- matrix(nrow=reps,ncol=size)
  
  for (i in 1:reps) {
    cmrep <- rrarefy(communityMatrix,sampleSize)
    for (j in 1:size) {
      rdt[i,j] <- d(cmrep,lev=levels[j],q=qs[j])
    }
  }
  
  meanD <- apply(rdt,MARGIN=2, FUN=mean)
  sterrD <- apply(rdt,MARGIN=2, FUN=sd)/sqrt(reps)
  
  diversity_table <- data.frame(row.names=paste(levels, qs, sep=""))
  
  diversity_table$mean <- meanD
  diversity_table$se <- sterrD
  #diversity_table$size <- sampleSize
  
  return (diversity_table)
}

# write rarefaction table file
writeRdiversityTable <- function(communityMatrix, outputRFTable, rareReps = 10, sampleSizesPoints = 11) {
  raremax <- min(rowSums(communityMatrix))
  #  rareReps <- 10
  #  sampleSizesPoints <- 11 # include 1
  
  sampleSizesSeq <- round(seq(1, raremax, length.out = sampleSizesPoints), digits = 0)
  
  rdiversityTableMean <- NULL
  rdiversityTableSe <- NULL
  for(sS in sampleSizesSeq){   
    if (is.null(rdiversityTableMean)) {
      rdiversityTableMean <- data.frame(row.names=paste(levels, qs, sep=""), check.names=FALSE)
      rdiversityTableMean <- rdiversityTable(communityMatrix, sS, rareReps)[,1]
      rdiversityTableSe <- data.frame(row.names=paste(levels, qs, sep=""), check.names=FALSE)
      rdiversityTableSe <- rdiversityTable(communityMatrix, sS, rareReps)[,2]
    } else {
      tmpRDT <- rdiversityTable(communityMatrix, sS, rareReps)[,1] 
      rdiversityTableMean <- cbind(rdiversityTableMean, tmpRDT)
      tmpRDT <- rdiversityTable(communityMatrix, sS, rareReps)[,2] 
      rdiversityTableSe <- cbind(rdiversityTableSe, tmpRDT)
    }  
    
    cat("Rarefaction take sample size ", sS, ", max = ", sampleSizesSeq[length(sampleSizesSeq)], "\n") 
  }
  colnames(rdiversityTableMean) <- paste("size.", sampleSizesSeq, sep="")
  colnames(rdiversityTableSe) <- paste("size.", sampleSizesSeq, sep="")
  
  write.csv(rdiversityTableMean, outputRFTable, row.names=FALSE, quote=FALSE)
  
  #LM.rdt <- lm(rdiversityTableMean[1,]~sampleSizesSeq, rdiversityTableMean[1,])
  
  cat("write rarefaction table file to ", outputRFTable, "\n")
}

