
# change config below
#sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
#setwd(sourcePath)
#workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
#matrixNames <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name   

if(!exists("tableFile")) stop("table file is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("isPlot")) isPlot = TRUE # by plot
if(!exists("otuThr")) otuThr = 97
if(!exists("levels")) levels = rep(c("gamma","alpha","beta"),3)
if(!exists("qs")) qs = rep(0:2,each=3)

n <- length(matrixNames) 

otusRowNames <- c("Reads", "OTUs", "Singleton", "Coupleton") # cannot get unique reads from CM
divRowNames <- c("$^0D_\\gamma$","$^0D_\\alpha$","$^0D_\\beta$","$^1D_\\gamma$","$^1D_\\alpha$","$^1D_\\beta$",
				"$^2D_\\gamma$","$^2D_\\alpha$","$^2D_\\beta$")
lotus = length(otusRowNames)

#source("Modules/init.R", local=TRUE)

######## summary by datasets including singletons #######

cat("\nTable: summary by datasets including singletons: isPlot =", isPlot, ", otuThr =", otuThr, "\n") 

# last column total
by.datasets <- data.frame(row.names=c(otusRowNames, divRowNames), stringsAsFactors=FALSE)

for (expId in 1:n) {	
  communityMatrix <- getCommunityMatrixT(expId, isPlot, FALSE)
  
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
  communityMatrix <- getCommunityMatrixT(expId, isPlot, rmSingleton)
  
  # TODO: 
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
               label = postfix(paste("tab:perSample", matrixNames[expId], sep=":"), isPlot, rmSingleton, sep=":"), 
               caption.placement = "top"), sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)
}
