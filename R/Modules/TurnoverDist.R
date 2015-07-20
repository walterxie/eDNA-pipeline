library(vegetarian)

TurnoverDist<-function(input.table){
  # Returns a distance matrix composed of pair-wise turnovers
  # Input: a community matrix
  # Depends on:
  #   vegetarian library 
   
  to.table<-matrix(0,nrow=nrow(input.table),ncol=nrow(input.table))
  for(i in 1:nrow(input.table)){
    for(j in 1:nrow(input.table)){
      to.table[i,j]<-turnover(input.table[c(i,j),])
    }
  }
	
  d <- as.dist(to.table)
  attr(d, "Labels") <- dimnames(input.table)[[1L]]
	
  return(d)
}