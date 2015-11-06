

if(!exists("tableFile")) stop("table file is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("isPlot")) isPlot = TRUE # by plot
if(!exists("otuThr")) otuThr = 97
if(!exists("levels")) levels = rep(c("gamma","alpha","beta"),3)
if(!exists("qs")) qs = rep(0:2,each=3)

n <- length(matrixNames) 

#source("Modules/init.R", local=TRUE)

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

  align.v <- rep("r", ncol(perSample) + 1)
  printXTable(perSample, caption = paste("Table of biodiversities per sample for", matrixNames[expId]),
              label = postfix(paste("tab:perSample", matrixNames[expId], sep=":"), isP, min2, sep=":"), 
              align = align.v, file=tableFile)
}
