
library(vegan)

if(!exists("levels")) 
	levels = rep(c("gamma","alpha","beta"),3)
if(!exists("qs")) 
	qs = rep(0:2,each=3)

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




