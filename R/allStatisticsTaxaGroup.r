

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

cat("\nTable: summary by datasets by taxa groups: rmSingleton =", rmSingleton, "isPlot =", isPlot, ", otuThr =", otuThr, "\n") 

# last column total
by.datasets <- data.frame(row.names=taxa.groups, stringsAsFactors=FALSE)

for (expId in 1:(n-1)) {	
  isP <- isPlot
  if (expId == n) {
    isP <- TRUE
  } 
  communityMatrix <- getCommunityMatrixT(expId, isP, FALSE) # always including singletons
  
  diversity_table <- diversity.df(communityMatrix)
  
  by.datasets[1,expId] <- sum(communityMatrix)
  by.datasets[2,expId] <- ncol(communityMatrix)
  by.datasets[3,expId] <- sum(colSums(communityMatrix)==1)
  by.datasets[4,expId] <- sum(colSums(communityMatrix)==2)
  for (i in 1:length(unlist(diversity_table))) {
    by.datasets[i+lotus,expId] <- unlist(diversity_table)[i]
  }
  colnames(by.datasets)[expId] <- matrixNames[expId]
}

by.datasets[1:4,"Total"] <- rowSums(by.datasets)[1:4]
#by.datasets[1:4,] <- round(by.datasets[1:4,], 0)
#by.datasets[5:nrow(by.datasets),] <- round(by.datasets[5:nrow(by.datasets),], 2)
by.datasets <- round(by.datasets, 2)
by.datasets <- format(by.datasets, big.mark=",", scientific=F)
for (i in 1:lotus)
  by.datasets[i,] <- gsub(".00", "", by.datasets[i,])

print(xtable(by.datasets, caption = paste("Table of sequence statistics for eDNA data sets", 
                                          paste(matrixNames, collapse = ", ")),
             label = "tab:biodiveDNA", caption.placement = "top"), 
             sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

######## summary by plots #######

cat("\nTable: summary by plots/subplot: isPlot =", isPlot, "rmSingleton =", rmSingleton, ", otuThr =", otuThr, "\n") 

for (expId in 1:n) {	
  isP <- isPlot
  min2 <- rmSingleton
  if (expId == n) {
    isP <- TRUE
    min2 <- FALSE
  } 
  communityMatrix <- getCommunityMatrixT(expId, isP, min2)
  source("Modules/Diversities.R")
  
  # reads
  readsPerSample <- abundancePerSample(communityMatrix, hasTotal=TRUE) 
  # OTUs
  otusPerSample <- richnessPerSample(communityMatrix, hasTotal=TRUE)
  # gamme1
  shanPerSample <- shannonPerSample(communityMatrix, digits=2)
  
  if (!any( tolower(rownames(readsPerSample)) == tolower(rownames(otusPerSample)) )) 
    stop("OTU names do not match !")
  if (!any( tolower(rownames(readsPerSample)[! rownames(readsPerSample) %in% "Total"]) 
            == tolower(rownames(shanPerSample)) )) 
    stop("OTU names do not match !")
  
  perSample <- data.frame(reads=readsPerSample, OTUs=otusPerSample, stringsAsFactors=FALSE)
  
  shanPerSample <- rbind(shanPerSample, c("NA"))
  perSample<- cbind(perSample, shanPerSample)
  colnames(perSample) <- c("reads","OTUs","Shannon")
  
#  if (sort=="descending") {
    perSample <- perSample[order(perSample[,"reads"], decreasing = TRUE), ]
#  } else if (sort=="ascending") {
#    perSample <- perSample[order(perSample$abundance), ]
#  }# sort=="no"
  
  perSample <- format(perSample, big.mark=",", scientific=F)

  print(xtable(perSample, caption = paste("Table of biodiversities per sample for", matrixNames[expId]),
               label = postfix(paste("tab:perSample", matrixNames[expId], sep=":"), isP, min2, sep=":"), 
               caption.placement = "top"), sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)
}
