# "-by-plot" trigger merge 2 subplots columns
library(gplots)
library(ggplot2)
library(grid)
library(RColorBrewer)

#figDir <- "figures"
#sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
#setwd(sourcePath)
#workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
#matrixNames <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name   
#levels = rep(c("gamma","alpha","beta"),3)
#qs = rep(0:2,each=3)
#rmSingleton <- FALSE 

if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("rmSingleton")) stop("rmSingleton flag is missing !")

otuThr = 97

subTitles <- c("(a)","(b)","(c)","(d)","(e)","(f)") # "(a)","(b)","(c)","(d)","(e)","(f)"
subTitles <- paste(subTitles, matrixNames, sep = " ")

# show the phyla that have at least percThr of maximum OTUs in each marker, and group all the rest as "Other".
percThr = 0.001

# http://stackoverflow.com/questions/22295253/force-bars-to-start-from-a-lower-value-than-0-in-ggplot-geom-bar-in-r
# defining the scale change
# scale_y_continuous(trans = mylog_trans(base=10, from=-2)) # starts from 1e-2
library(scales)
mylog_trans <- function(base=exp(1), from=0) {
  trans <- function(x) log(x, base)-from
  inv <- function(x) base^(x+from)
  trans_new("mylog", trans, inv, log_breaks(base=base), domain = c(base^from, Inf))
}

taxaAssgFile <- paste(workingPath, "data/taxonomy97phyla.txt", sep="") 

col.names <- readLines(taxaAssgFile, n=1)
col.names <- unlist(strsplit(col.names, split="\t")) # 1st is "Phylumn" 
# readLines seems to change file pointer
taxaAssg <- read.csv(taxaAssgFile, header=T, row.names=1, sep = "\t", col.names=col.names, check.names=F)  # last column is taxagroup 

nc <- ncol(taxaAssg) - 1
otuSum <- colSums(taxaAssg[,-ncol(taxaAssg)])
nr <- nrow(taxaAssg)

print(paste("Total phyla = ", nr, sep=""))  
for (co in 1:nc) {
	print(paste(names(taxaAssg)[co], "'s phyla = ", length(which(taxaAssg[,co] != 0)), ", total OTUs = ", sum(taxaAssg[,co]), sep=""))
}

# rm singleton OTU from taxa
if (rmSingleton) {
	allTaxaSingleton <- read.table(paste(workingPath, "data/taxonomy97singleton.txt", sep=""), sep="\t", row.names=1, header=T, check.names=F)
	print("Find singleton OTU:")
	print(colSums(allTaxaSingleton))
	
	matched <- match(rownames(allTaxaSingleton), rownames(taxaAssg))

	for (r in 1:nrow(taxaAssg)) {
		# has singleton OTU
		if (is.element(r, matched)) 
			taxaAssg[r,-7] <- taxaAssg[r,-7] - allTaxaSingleton[which(matched == r),]
	}
	
	# rm 0 rows
	taxaAssg <- taxaAssg[-which(rowSums(taxaAssg[,-7])==0),]
	
	nc <- ncol(taxaAssg) - 1
	otuSum <- colSums(taxaAssg[,-ncol(taxaAssg)])
	nr <- nrow(taxaAssg)
	
	print(paste("Total left phyla = ", nr, sep="")) 
	for (co in 1:nc) {
		print(paste(names(taxaAssg)[co], "'s phyla = ", length(which(taxaAssg[,co] != 0)), ", total OTUs = ", sum(taxaAssg[,co]), sep=""))
	} 
}

#colSums(taxaAssgPer)==1
taxaAssgPer <- prop.table(as.matrix(taxaAssg[,-ncol(taxaAssg)]), 2)

taxaAssgOTUs <- data.frame(taxaAssg, check.names=F)
taxaAssgOTUs <- taxaAssgOTUs[apply(taxaAssgPer, 1, max) >= percThr,]
# Only show the phyla that are at least 0.1% in each, and group all the rest as "Other".
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

print(paste(names(taxaAssgOTUs)[1], "'s selected phyla = ", length(which(taxaAssgOTUs[,1] != 0)), ", identified OTUs = ", sum(taxaAssgOTUseDNA$OTUs), sep=""))

for (co in 2:nc) {
	tmpeDNA <- data.frame(phylumn= row.names(taxaAssgOTUs), OTUs= taxaAssgOTUs[,co], eDNA= names(taxaAssgOTUs)[co], taxagroup = taxaAssgOTUs[,ncol(taxaAssgOTUs)], stringsAsFactors = FALSE)
	if (sum(other)>0) {
		otherOTUs <- data.frame(phylumn= "Others", OTUs= other[co] * otuSum[co], eDNA= names(taxaAssgOTUs)[co], taxagroup = "Others", stringsAsFactors = FALSE)
		tmpeDNA <- rbind(tmpeDNA, otherOTUs)
	}
	print(paste(names(taxaAssgOTUs)[co], "'s selected phyla = ", length(which(taxaAssgOTUs[,co] != 0)), ", total OTUs = ", sum(tmpeDNA$OTUs), sep=""))

	taxaAssgOTUseDNA <- rbind(taxaAssgOTUseDNA, tmpeDNA)
} 
rownames(taxaAssgOTUseDNA) <- NULL

print(paste("Selected phyla = ", nrow(taxaAssgPer), sep=""))  

taxaAssgOTUseDNA$eDNA = factor(taxaAssgOTUseDNA$eDNA, names(taxaAssgOTUs)[-ncol(taxaAssgOTUs)])
taxaAssgOTUseDNA$phylumn = factor(taxaAssgOTUseDNA$phylumn, rev(tmpeDNA$phylumn))
taxaAssgOTUseDNA$taxagroup = factor(taxaAssgOTUseDNA$taxagroup, unique(taxaAssgOTUseDNA$taxagroup))

###### phyla OTUs bar #####

if (percThr == 0) {
	xlab <- xlab(paste(nrow(taxaAssgPer), " phyla + high-level taxa in total", sep=""))
} else {
	xlab <- xlab(paste(nrow(taxaAssgPer), " phyla + high-level taxa (OTUs >= ", percThr*100, "% of maximum)", sep=""))
}

pdf(paste(workingPath, figDir, "/all-taxa-phyla-OTUs-", otuThr, ".pdf", sep = ""), width=12, height=11)	

print( ggplot(taxaAssgOTUseDNA, aes(x = phylumn, y = OTUs, fill = taxagroup)) + geom_bar(stat = "identity", position = "identity") + 
	scale_y_continuous(trans = mylog_trans(base=10, from=-0.3), labels = function (x) floor(x), 
		breaks=c(0.1, 1, 10, 100, 1000, 10000), expand = c(0, 0)) +
     coord_flip() + theme_bw() + facet_grid( ~ eDNA) + ylab("OTUs") + xlab + 
	 theme(legend.position=c(0.4, 1.04), legend.direction="horizontal", plot.margin=unit(c(1.2,0.5,0.2,0.2), "cm"), panel.margin = unit(0.8, "lines")) + 
	 scale_fill_brewer(name = "Taxa Group : ",palette = "Paired") )
	 
invisible(dev.off()) 

