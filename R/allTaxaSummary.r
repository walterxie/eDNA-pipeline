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
# plot bar chart either for reads or OTUs, determined by taxaAssg created by getTaxaAssg(y_string=?)
source("Modules/TaxonomyAssignment.R")

y_string="reads" # or OTUs

######## eDNA by plot/subplot #########
taxaAssgTotal <- NULL
for (expId in 1:(n-1)) {
  taxaAssgReads <- getTaxaAssgReads(expId, isPlot, rmSingleton, rankLevel, groupLevel, taxa.group) 
  
  if (nrow(taxaAssgReads) < 1) {
    cat("\nSkip", taxa.group, "subset taxonomic summary from", matrixNames[expId], ".\n") 
    next
  }
  
  taxaAssg <- getTaxaAssg(taxaAssgReads, rankLevel, groupLevel, y_string) 
  
  cat(matrixNames[expId], ":", rankLevel, "=", nrow(taxaAssg), ",", groupLevel, "=", length(unique(taxaAssg[,groupLevel])), ".\n")
  
  if (expId==1) {
    taxaAssgTotal <- taxaAssg[,c(rankLevel,total_string,groupLevel)]
  } else {
    taxaAssgTotal <- merge(taxaAssgTotal, taxaAssg[,c(rankLevel,total_string,groupLevel)], 
                           by = c(rankLevel, groupLevel), all = TRUE)
  }
  colnames(taxaAssgTotal) <- gsub(total_string, matrixNames[expId], colnames(taxaAssgTotal))
  
  legend_nrow=1
  tA <- taxonomyAssignment(taxaAssg, rankLevel, groupLevel, percThr, legend_nrow, y_string, plotTotal=FALSE)
  pdfWidth = 0.1 + tA$maxLabelLen / 10 + (tA$ncol-2) * 1.5
  pdfHeight = 1 + legend_nrow * 1 + tA$nrow * 0.12
  
  fname <- paste(rankLevel, matrixNames[expId], postfix(taxa.group, TRUE, rmSingleton, sep="-"), y_string, sep = "-")
  pdf(file.path(workingPath, figDir, paste(fname, "pdf", sep = ".")), width=pdfWidth, height=pdfHeight)	
  print(tA$plot)
  invisible(dev.off()) 
}

###### total by gene #####
taxaAssgTotal[is.na(taxaAssgTotal)] <- 0
taxaAssgTotal <- taxaAssgTotal[order(taxaAssgTotal[, groupLevel]),]

cat(total_string, ":", rankLevel, "=", nrow(taxaAssgTotal), ",", groupLevel, "=", length(unique(taxaAssgTotal[,groupLevel])), ".\n")  

# add Total column
taxaAssgTotal[,total_string] <- rowSums(taxaAssgTotal[,-1:-2]) 

# move groupLevel to the last col
taxaAssgTotal <- taxaAssgTotal[,c(1,3:ncol(taxaAssgTotal),2)]

legend_nrow=1
tA <- taxonomyAssignment(taxaAssgTotal, rankLevel, groupLevel, percThr, legend_nrow, y_string, plotTotal=FALSE)
pdfWidth = 0.1 + tA$maxLabelLen / 10 + (tA$ncol-2) * 1.5
pdfHeight = 1 + legend_nrow * 1 + tA$nrow * 0.12

fname <- paste(rankLevel, tolower(total_string), postfix(taxa.group, TRUE, rmSingleton, sep="-"), y_string, sep = "-")
pdf(file.path(workingPath, figDir, paste(fname, "pdf", sep = ".")), width=pdfWidth, height=pdfHeight)	
print(tA$plot)
invisible(dev.off()) 
