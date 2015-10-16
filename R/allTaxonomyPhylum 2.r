library(gplots)
library(ggplot2)
library(grid)
library(RColorBrewer)


if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97

n <- length(matrixNames) 

# no domain of "prokaryote", Bacteria's parent is "cellular organisms", and same as Archaea
rankLevel="phylum" 
groupLevel="kingdom" # gives colour, and must higher than rankLevel
y_string="reads" # or OTUs

# show the phyla that have at least percThr of maximum OTUs in each marker, and group all the rest as "Other".
percThr = 0

source("Modules/init.r")
source("Modules/TaxonomyAssignment.R")

######## eDNA #########
for (expId in 1:(n-1)) {
  taxaAssgReads <- getTaxaAssgReads(expId, rmSingleton, rankLevel, groupLevel) 
    
  taxaAssg <- getTaxaAssg(taxaAssgReads, rankLevel, y_string) 
  
  legend_nrow=1
  tA <- taxonomyAssignment(taxaAssg, rankLevel, groupLevel, percThr, legend_nrow, y_string, plotTotal=FALSE)
  pdfWidth = 0.1 + tA$maxLabelLen / 10 + (tA$ncol-2) * 1.5
  pdfHeight = 1 + legend_nrow * 1 + tA$nrow * 0.12
  
  pdf(paste(workingPath, "figures/all-taxa-", rankLevel, "-", y_string, ".pdf", sep = ""), width=pdfWidth, height=pdfHeight)	
  print(tA$plot)
  invisible(dev.off()) 
}

col.names <- readLines(taxaAssgFile, n=1)
col.names <- unlist(strsplit(col.names, split="\t")) # 1st is "Phylumn" 
# readLines seems to change file pointer
taxaAssg <- read.csv(taxaAssgFile, header=T, row.names=1, sep = "\t", col.names=col.names, check.names=F)  # last column is taxagroup 
taxaAssg[is.na(taxaAssg)] <- 0

nc <- ncol(taxaAssg) - 1
otuSum <- colSums(taxaAssg[,-ncol(taxaAssg)])
nr <- nrow(taxaAssg)

print(paste("Total phyla = ", nr, sep=""))  

taxaAssgPer <- prop.table(as.matrix(taxaAssg[,-ncol(taxaAssg)]), 2)
taxaAssgOTUs <- data.frame(taxaAssg, check.names=F)
taxaAssgOTUs <- taxaAssgOTUs[apply(taxaAssgPer, 1, max) >= percThr,]
# Only show the phyla that are at least 1% in each, and group all the rest as "Other".
other <- colSums(taxaAssgPer[apply(taxaAssgPer, 1, max) < percThr,])
taxaAssgPer <- taxaAssgPer[apply(taxaAssgPer, 1, max) >= percThr,]

if (sum(other)>0) {
	taxaAssgPer <- rbind(taxaAssgPer, other)
	row.names(taxaAssgPer)[nrow(taxaAssgPer)] <- "Others"  
}

taxaAssgOTUseDNA <- data.frame(phylumn= row.names(taxaAssgOTUs), OTUs= taxaAssgOTUs[,1], eDNA= names(taxaAssgOTUs)[1], taxagroup = taxaAssgOTUs[,ncol(taxaAssgOTUs)], stringsAsFactors = FALSE)
if (sum(other)>0) {
	otherOTUs <- data.frame(phylumn= "Others", OTUs= other[1] * otuSum[1], eDNA= names(taxaAssgOTUs)[1], taxagroup = "Others", stringsAsFactors = FALSE)
	taxaAssgOTUseDNA <- rbind(taxaAssgOTUseDNA, otherOTUs)
}

print(paste(names(taxaAssgOTUs)[1], "'s phyla = ", length(which(taxaAssgOTUs[,1] != 0)), ", identified OTUs = ", sum(taxaAssgOTUseDNA$OTUs), sep=""))

for (co in 2:nc) {
	tmpeDNA <- data.frame(phylumn= row.names(taxaAssgOTUs), OTUs= taxaAssgOTUs[,co], eDNA= names(taxaAssgOTUs)[co], taxagroup = taxaAssgOTUs[,ncol(taxaAssgOTUs)], stringsAsFactors = FALSE)
	if (sum(other)>0) {
		otherOTUs <- data.frame(phylumn= "Others", OTUs= other[co] * otuSum[co], eDNA= names(taxaAssgOTUs)[co], taxagroup = "Others", stringsAsFactors = FALSE)
		tmpeDNA <- rbind(tmpeDNA, otherOTUs)
	}
	print(paste(names(taxaAssgOTUs)[co], "'s phyla = ", length(which(taxaAssgOTUs[,co] != 0)), ", total OTUs = ", sum(tmpeDNA$OTUs), sep=""))

	taxaAssgOTUseDNA <- rbind(taxaAssgOTUseDNA, tmpeDNA)
} 
rownames(taxaAssgOTUseDNA) <- NULL

print(paste("Selected phyla = ", nrow(taxaAssgPer), sep=""))  

taxaAssgOTUseDNA$eDNA = factor(taxaAssgOTUseDNA$eDNA, names(taxaAssgOTUs)[-ncol(taxaAssgOTUs)])
taxaAssgOTUseDNA$phylumn = factor(taxaAssgOTUseDNA$phylumn, rev(tmpeDNA$phylumn))
taxaAssgOTUseDNA$taxagroup = factor(taxaAssgOTUseDNA$taxagroup, unique(taxaAssgOTUseDNA$taxagroup))

###### phyla OTUs bar #####

pdf(paste(workingPath, "figures/all-taxa-phyla-OTUs", postfix, ".pdf", sep = ""), width=12, height=11)	

print( ggplot(taxaAssgOTUseDNA, aes(x = phylumn, y = OTUs, fill = taxagroup)) + geom_bar(stat = "identity", position = "identity") + 
     coord_flip() + scale_y_log10(breaks=c(1, 10, 100, 1000, 10000), expand = c(0, 0)) +
	 theme_bw() + facet_grid( ~ eDNA) + ylab("OTUs") + xlab(paste(nrow(taxaAssgPer), " phyla + equivalent-level taxa (OTUs >= ", percThr*100, "% of maximum)", sep="")) + #expression("phyla (maximum OTUs" >= " 0.1%))"
	 theme(legend.position=c(0.4, 1.04), legend.direction="horizontal", plot.margin=unit(c(1.2,0.5,0.2,0.2), "cm"), panel.margin = unit(0.8, "lines")) + 
	 scale_fill_brewer(name = "Taxa Group : ",palette = "Paired") )

invisible(dev.off()) 

