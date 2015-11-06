

if(!exists("tableFile")) stop("table file is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")
if(!exists("taxa.groups")) stop("taxa groups are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("isPlot")) isPlot = TRUE # by plot
if(!exists("otuThr")) otuThr = 97
if(!exists("levels")) levels = rep(c("gamma","alpha","beta"),3)
if(!exists("qs")) qs = rep(0:2,each=3)


n <- length(matrixNames) 

#source("Modules/init.R", local=TRUE)

######## summary by datasets including singletons #######
otusRowNames <- c("Reads", "OTUs", "Singleton", "Coupleton") # cannot get unique reads from CM

cat("\nTable: summary by datasets by taxa groups: rmSingleton =", rmSingleton, ", otuThr =", otuThr, "\n") 

# last column total
df.otus <- data.frame(row.names=taxa.groups, stringsAsFactors=FALSE)
df.reads <- data.frame(row.names=taxa.groups, stringsAsFactors=FALSE)

for (expId in 1:(n-1)) {	
  n.taxa <- 0
  for (taxa.group in taxa.groups) {
    cat("Summary by", taxa.group, "subset from", matrixNames[expId], ". \n")
    
    # always by plot here to save time
    communityMatrix <- getCommunityMatrixT(expId, TRUE, rmSingleton, taxa.group)
    
    if (is.null(communityMatrix)) {
      cat("Warning:", taxa.group, "subset from", matrixNames[expId], "has no OTU.\n") 
    } 
    
    n.taxa <- n.taxa + 1
    if (ncol(communityMatrix) > 0) {
      df.otus[n.taxa,expId] <- ncol(communityMatrix)
      df.reads[n.taxa,expId] <- sum(communityMatrix)
    } else {
      df.otus[n.taxa,expId] <- 0
      df.reads[n.taxa,expId] <- 0
    }
    colnames(df.otus)[expId] <- matrixNames[expId]
    colnames(df.reads)[expId] <- matrixNames[expId]
  }
}
# add total
df.otus[,"Total"] <- rowSums(df.otus)
df.reads[,"Total"] <- rowSums(df.reads)

df.otus <- format(df.otus, big.mark=",", scientific=F)
df.reads <- format(df.reads, big.mark=",", scientific=F)

align.v <- rep("r", ncol(df.otus) + 1)
printXTable(df.otus, caption = paste("Table of OTUs statistics for eDNA data sets by taxonomic groups"), 
            align = align.v, label = paste("tab:otus","all", sep = ":"), file=tableFile)

align.v <- rep("r", ncol(df.reads) + 1)
printXTable(df.reads, caption = paste("Table of reads statistics for eDNA data sets by taxonomic groups"), 
            align = align.v, label = paste("tab:reads", "all", sep = ":"), file=tableFile)


