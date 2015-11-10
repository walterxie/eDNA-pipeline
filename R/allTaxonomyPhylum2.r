library(gplots)
library(ggplot2)
library(grid)
library(RColorBrewer)


if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("isPlot")) isPlot = FALSE # by subplot
if(!exists("taxa.group")) taxa.group="all"
if(!exists("rankLevel")) rankLevel="phylum"
if(!exists("groupLevel")) groupLevel="kingdom"

n <- length(matrixNames) 

total_string = "Total"
# no domain of "prokaryote", Bacteria's parent is "cellular organisms", and same as Archaea
#rankLevel="phylum" 
#groupLevel="kingdom" # gives colour, and must higher than rankLevel
# show the phyla that have at least percThr of maximum OTUs in each marker, and group all the rest as "Other".
percThr = 0

source("Modules/init.r")
# plot OTUs and reads in one graph
source("Modules/TaxonomyAssignment2.R")

######## eDNA by plot/subplot #########
taxaAssgTotal <- NULL
for (expId in 1:(n-1)) {
  taxaAssgReads <- getTaxaAssgReads(expId, isPlot, rmSingleton, rankLevel, groupLevel, taxa.group) 
  
  if (nrow(taxaAssgReads) < 1) {
    cat("\nSkip", taxa.group, "subset taxonomic summary from", matrixNames[expId], ".\n") 
    next
  }
  
  taxaAssg <- getTaxaAssg2(taxaAssgReads, rankLevel, groupLevel) 
  
  cat(matrixNames[expId], ":", rankLevel, "=", nrow(taxaAssg), ",", groupLevel, "=", length(unique(taxaAssg[,groupLevel])), ".\n")
  
  if (expId==1) {
    taxaAssgTotal <- taxaAssg[,c(rankLevel,total_string,groupLevel,"var")]
  } else {
    taxaAssgTotal <- merge(taxaAssgTotal, taxaAssg[,c(rankLevel,total_string,groupLevel,"var")], 
                           by = c(rankLevel, groupLevel,"var"), all = TRUE)
  }
  colnames(taxaAssgTotal) <- gsub(total_string, matrixNames[expId], colnames(taxaAssgTotal))
 
  legend_nrow=1
  tA <- taxonomyAssignment2(taxaAssg, rankLevel, groupLevel, percThr, legend_nrow, plotTotal=FALSE)
  pdfWidth = 0.1 + tA$maxLabelLen / 10 + (tA$ncol-2) * 1.5
  pdfHeight = 1 + legend_nrow * 1 + tA$nrow * 0.8
  
  fname <- paste(rankLevel, matrixNames[expId], postfix(taxa.group, TRUE, rmSingleton, sep="-"), sep = "-")
  pdf(paste(workingPath, figDir, "/", fname, ".pdf", sep = ""), width=pdfWidth, height=pdfHeight)	
  print(tA$plot)
  invisible(dev.off()) 
}

###### total by gene #####
taxaAssgTotal[is.na(taxaAssgTotal)] <- 0

cat(total_string, ":", rankLevel, "=", nrow(taxaAssgTotal), ",", groupLevel, "=", length(unique(taxaAssgTotal[,groupLevel])), ".\n")  

# add Total column, without phylum,kingdom,var
taxaAssgTotal[,total_string] <- rowSums(taxaAssgTotal[,-1:-3]) 

# move groupLevel and var to the last col
taxaAssgTotal <- taxaAssgTotal[,c(1,4:ncol(taxaAssgTotal),2,3)]

legend_nrow=1
tA <- taxonomyAssignment2(taxaAssgTotal, rankLevel, groupLevel, percThr, legend_nrow, plotTotal=FALSE)
pdfWidth = 0.1 + tA$maxLabelLen / 10 + (tA$ncol-2) * 1.5
pdfHeight = 1 + legend_nrow * 1 + tA$nrow * 0.09

fname <- paste(rankLevel, tolower(total_string), postfix(taxa.group, TRUE, rmSingleton, sep="-"), sep = "-")
pdf(paste(workingPath, figDir, "/", fname, ".pdf", sep = ""), width=pdfWidth, height=pdfHeight)	
print(tA$plot)
invisible(dev.off()) 
