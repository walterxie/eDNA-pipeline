
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


