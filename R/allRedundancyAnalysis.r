# A Caution Regarding Rules of Thumb for Variance Inflation Factors
# http://link.springer.com/article/10.1007/s11135-006-9018-6

library(vegan)
library(ggplot2)
library(grid)
library(VIF)
library(xtable)
library(scales)

if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("tableFile")) stop("table file is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("diss.fun")) diss.fun="beta1-1"
if(!exists("taxa.group")) taxa.group="assigned"
if(!exists("isPlot")) isPlot = FALSE # by subplot

n <- length(matrixNames) 

source("Modules/RedundancyAnalysis.R", local=TRUE)

# Prepare environmental data for analysis -------------------------------------
env <- getSampleMetaData(isPlot)
# remove columns: Plot, Elev.group
env <- env[,-c(1,9)]
env[,"ForestType"] <- gsub(":.*", "", env[,"ForestType"], ignore.case = T)
env[env=="x"] <- NA
# CM30b51, CM30b58 are removed
env <- na.omit(env)
# remove discrete values
env <- env[,-c(4,7)]

# Inspect data
#plot(env, gap = 0, panel = panel.smooth)
# Log transform chem variables according inspection
env[,6:13] <- log(env[,6:13])
# log transformation for 0
env[env=="-Inf"] <- 0

# Inspect data again; Note that certain variables are highly correlated (e.g. EC/Organic.C/Total.N, 
# Elevation/Temperature). Also, two datasets are included - pilot vs full LBI study).
pdf(paste(workingPath, figDir, "/environmental-data.pdf", sep = ""), width=12, height=12)
plot(env, gap = 0, panel = panel.smooth, upper.panel=NULL) 
invisible(dev.off()) 

colnames(env)

#rdaReducedList <- list()
#rdaForwardList <- list()
#rdaBackwardList <- list()
#rdaDataList <- list()
#count = 1

minOTUs <- 200

for (expId in 1:(n-1)) {	
  cat("RDA Analysis:", taxa.group, "subset from", matrixNames[expId], 
      ", isPlot =", isPlot, ", rmSingleton =", rmSingleton, ", minOTUs =", minOTUs, ". \n")

  communityMatrix <- getCommunityMatrixT(expId, isPlot, rmSingleton, taxa.group, minRow=minOTUs)
  
  if (is.null(communityMatrix)) {
    cat("\nSkip", taxa.group, "subset from", matrixNames[expId], ".\n") 
    next
  }
  
  # merge to match rownames
  ncol.cm <- ncol(communityMatrix)
  cm.env <- merge(communityMatrix, env, by = "row.names")
  # move 1st col Row.names to row.names
  rownames(cm.env) <- cm.env[,"Row.names"]
  cm.env <- cm.env[,-1]
  
  communityMatrix <- cm.env[,1:ncol.cm]
  cm.env <- as.data.frame(data.matrix(cm.env[,-(1:ncol.cm)]))
  
  # remove 0 row/column after merge
  communityMatrix <- prepCommunityMatrix(communityMatrix)
  
  cat("RDA of", taxa.group, "subset from", matrixNames[expId], "having", ncol(communityMatrix), "OTUs", 
      nrow(communityMatrix), "samples. \n") 
  
  rda <- proceedRDA(cm.or.dist=communityMatrix, env=cm.env, tableFile, verbose)
    
  fname <- paste("rda", matrixNames[expId], postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-")
  pdf(file.path(workingPath, figDir, paste(fname, "pdf", sep = ".")), width=9, height=3)
  attach(mtcars)
  par(mfrow=c(1,3))
  par(mar=c(4,4,3,2), cex=0.6)
  
  rownames(rda$rda_reduced$CCA$wa) <- substrRight(rownames(rda$rda_reduced$CCA$wa), 4)
  rownames(rda$rda_forward$CCA$wa) <- substrRight(rownames(rda$rda_forward$CCA$wa), 4)
  rownames(rda$rda_backward$CCA$wa) <- substrRight(rownames(rda$rda_backward$CCA$wa), 4)
  
  plot(rda$rda_reduced, display = c("wa", "bp"), main="Reduced (VIF >= 10)")
  plot(rda$rda_forward, display = c("wa", "bp"), ylab="", main="Forward")
  plot(rda$rda_backward, display = c("wa", "bp"), ylab="", main="Backward")
  
  invisible(dev.off())
} # END for expId


