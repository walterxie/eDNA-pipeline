# "-by-plot" trigger merge 2 subplots columns
library(gplots)
library(ggplot2)
library(grid)
library(RColorBrewer)

workingPath <- "~/Subversion/compevol/research/NZGenomicObservatory/Metabarcoding/"
matrixNames <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # "16S", "18S", "trnL", "ITS", "COI", "COI-spun"
subTitles <- c("(a)","(b)","(c)","(d)","(e)","(f)") # "(a)","(b)","(c)","(d)","(e)","(f)"
subTitles <- paste(subTitles, matrixNames, sep = " ")

otuThr = 97
percThr = 0.0001

taxaAssgFile <- paste(workingPath, "data/taxonomy97phyla.txt", sep="") 

col.names <- readLines(taxaAssgFile, n=1)
col.names <- unlist(strsplit(col.names, split="\t")) # 1st is "Phylumn" 
# readLines seems to change file pointer
taxaAssg <- read.csv(taxaAssgFile, header=T, row.names=1, sep = "\t", col.names=col.names, check.names=F)  # last column is taxagroup 

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

print(paste(names(taxaAssgOTUs)[1], "'s phyla = ", length(which(taxaAssgOTUs[,1] != 0)), ", total OTUs = ", sum(taxaAssgOTUseDNA$OTUs), sep=""))

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

pdf(paste(workingPath, "figures/all-taxa-phyla-OTUs-", otuThr, ".pdf", sep = ""), width=12, height=11)	

print( ggplot(taxaAssgOTUseDNA, aes(x = phylumn, y = OTUs, fill = taxagroup)) + geom_bar(stat = "identity", position = "identity") + 
     coord_flip() + scale_y_log10(breaks=c(1, 10, 100, 1000, 10000), expand = c(0, 0)) +
	 theme_bw() + facet_grid( ~ eDNA) + ylab("OTUs") + xlab(paste("phyla + high-level taxa (", nrow(taxaAssgPer), ")", sep="")) + #expression("phyla (maximum OTUs" >= " 0.1%))"
	 theme(legend.position=c(0.4, 1.05), legend.direction="horizontal", plot.margin=unit(c(1.2,0.5,0.2,0.2), "cm"), panel.margin = unit(0.8, "lines")) + 
	 scale_fill_brewer(name = "Taxa Group : ",palette = "Paired") )

invisible(dev.off()) 

