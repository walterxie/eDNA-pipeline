# Create a figure: the number of reads assigned to a pre-defined rank level and also high-level taxa, and coloured by another given rank level.
# Input: 1) taxonomy table from MEGAN, 1st column is OTUs, the last column is taxonomy, the middles are samples (plots/subplots);
# Input: 2) community matrix, but it is the transverse matrix of community matrix in vegan package;
# Input: 3) rankLevel, the number of reads assigned to;
# Input: 4) groupLevel, for colouring, and must higher than rankLevel;
# Output: 1) bar chart.
# 
# Author: Walter Xie
# Accessed on 16 June 2015

library(gplots)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(tidyr)

# mv row in data frame
mvRow <- function(taxaAssg, from, to) {
	if (from < to) 
		stop("from < to Not implemented !")
		
	if (from == to) {
		return(taxaAssg)
	} else if (to == 1) {
		return(taxaAssg[c(from,1:(from-1),(from+1):nrow(taxaAssg)),])
	} else if (from == nrow(taxaAssg)) {
		return(taxaAssg[c(1:(to-1),from,(to):(from-1)),])
	} else {
		return(taxaAssg[c(1:(to-1),from,(to):(from-1),(from+1):nrow(taxaAssg)),])
	}
}

# change config below
workingPath <- "~/svn/compevol/research/NZGenomicObservatory/MiSeq/DOCt1/"
matrixNames <-  c("16S", "18S", "26S", "ITS", "FolCO1", "ShCO1") # only for cm file name and folder name   
taxaFiles <- c("DOCx9_16S_R1_reads_OTUs_num_nt_taxatable.txt", "DOCx9_18S_PEARcontigs_ME1.0_OTUs_num_nt_taxatable.txt", 
				"DOCx9_26S_R1_reads_OTUs_num_nt_taxatable.txt", "DOCx9_ITS_PEARcontigs_ME1.0_OTUs_num_nt_taxatable.txt", 
				"DOCx9_FolCO1_R1_reads_OTUs_num_nt_taxatable.txt", "DOCx9_ShCO1_PEARcontigs_ME1.0_OTUs_num_nt_taxatable.txt")
# unclassified category
unclassTaxa <- c("Not assigned", "No hits", "cellular organisms", "root")

# how taxa to display, choose either maxNTaxa or percThr
#maxNTaxa=30 
percThr = 0.001

# MEGAN mapping file column order
# no domain of "prokaryote", Bacteria's parent is "cellular organisms", and same as Archaea
ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus")
rankLevel="order" 
groupLevel="phylum" # gives colour, and must higher than rankLevel

myPalette <- c("#CCCCCC", "#999999", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", 
			   "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#E69F00", "#FF6666", "#56B4E9", "#009E73", "#F0E442", "#9933FF",   
				"#0072B2", "#D55E00", "#CC79A7", "#E6AB02", "#33ff99", "#fb9a99", "#e31a1c", "#fdbf6f", "#cab2d6", "#e41a1c", 
				"#377eb8", "#ff9933", "#4daf4a", "#984ea3", "#cc0033", "#a65628", "#9999FF", "#f781bf", "#a6cee3", "#1f78b4", 
				"#b2df8a", "#a6cee3", "#ffff33", "#006600")


for (expId in 1:length(matrixNames)) {
    print(paste("Load data for", matrixNames[expId], " : "))
    
	##### load data #####
	inputCM <- paste(workingPath, "data/", matrixNames[expId], ".txt", sep="") 
	communityMatrix <- read.table(inputCM, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
	communityMatrix <- communityMatrix[order(rownames(communityMatrix)),]
	communityMatrix <- communityMatrix[,order(colnames(communityMatrix))]
	# U115 = U155 (correct name) and OJ95 = CJ95 (correct name). AL132 is the "difficult" site extracting DNA 
	colnames(communityMatrix) <- gsub("OJ95", "CJ95", colnames(communityMatrix), ignore.case = T)
	colnames(communityMatrix) <- gsub("U115", "U155", colnames(communityMatrix), ignore.case = T)
	
	inputTaxa <- paste(workingPath, "Taxonomy_tables/", taxaFiles[expId], sep="")
	taxaPaths <- read.table(inputTaxa, header=TRUE, row.names=1, sep = "\t", stringsAsFactors=FALSE)  
	taxaPaths <- taxaPaths[order(rownames(taxaPaths)),]

	rownames(communityMatrix) <- gsub(";size.*;", "", rownames(communityMatrix))
	rownames(taxaPaths) <- gsub(";size.*;", "", rownames(taxaPaths))

	if ( all( tolower(rownames(communityMatrix)) != tolower(rownames(taxaPaths)) ) ) 
		stop( paste("OTU names in", matrixNames[expId], "community matrix file and taxa path file are not matched !") )

	if ( ! rankLevel %in% colnames(taxaPaths) ) 
		stop( paste("Column name", rankLevel, "not exist in taxa path file for ", matrixNames[expId]) )
	if (! groupLevel %in% colnames(taxaPaths) ) 
		stop( paste("Column name", groupLevel, "not exist in taxa path file for ", matrixNames[expId]) )

	# define unclassified
	taxaPaths[taxaPaths==unclassTaxa[1] | taxaPaths==unclassTaxa[2] | taxaPaths==unclassTaxa[3] | is.na(taxaPaths)] <- "unclassified"
	
	###### taxa assignment by reads #####
	# create taxa table based on the number of reads
	taxaAssgReads <- cbind(communityMatrix, taxaPaths[,rankLevel])
	colnames(taxaAssgReads)[ncol(taxaAssgReads)] <- rankLevel

	taxaAssg <- aggregate(as.formula(paste(".", rankLevel, sep=" ~ ")), data= taxaAssgReads, FUN=sum)
	taxaAssg$Total <- rowSums(taxaAssg[,-1]) 

	# get col index for taxaAssg: 1 is rankLevel, 2:colTotal are number of reads by plots and total, colTotal+1:ncol are taxa breaks
	colRanks <- match(ranks, colnames(taxaPaths))
	colRankLevel <- which(colnames(taxaPaths)==rankLevel)
	colGroupLevel <- which(colnames(taxaPaths)==groupLevel)
	
	if (colRankLevel <= colGroupLevel || colGroupLevel <= colRanks[1]) 
		stop( paste("groupLevel (", groupLevel, ") must be higher than rankLevel (", rankLevel, 
			", and lower than colRanks[1] (", colRanks[1], ") !") )

	# add all higher taxa groups than rankLevel	
	for ( taxaCol in sort(colRanks[colRanks<colRankLevel], decreasing = T) ) {
		matchTaxa <- match(taxaAssg[,1], taxaPaths[,rankLevel])	
		taxaAssg[,colnames(taxaPaths)[taxaCol]] <- taxaPaths[matchTaxa,taxaCol]
	}
	
	#####  sort by each level of taxa group #####
	colTotal <- which(colnames(taxaAssg)=="Total")
	colTaxaGroupLevel <- which(colnames(taxaAssg)==groupLevel)
	
	if ((colTotal+1) >= ncol(taxaAssg)) 
		stop( paste("Taxa groups should have more than 2 columns : ", ncol(taxaAssg)-colTotal, " !") )

	# taxaAssg[order(taxaAssg[,15],taxaAssg[,14],taxaAssg[,13],taxaAssg[,12],-taxaAssg$Total),]
	# taxaAssg <- taxaAssg[order(taxaAssg$Total, decreasing = T),]
	taxaAssg <- eval(parse( text=paste("taxaAssg[order(taxaAssg[,", 
					paste(ncol(taxaAssg):(colTotal+1), collapse="],taxaAssg[,"), "],-taxaAssg$Total),]", sep="") ))
	rownames(taxaAssg) <- 1:nrow(taxaAssg)
	
	# 2nd last col taxagroup for colouring
	taxaAssg$taxagroup <- taxaAssg[,colTaxaGroupLevel]
	# last col space to show taxa label as rank hierarchy
	taxaAssg$space <- rep("", nrow(taxaAssg))
	
	# take "unclassified" row out temporarily for easy coding
	unclId <- which(taxaAssg[,1]=="unclassified")
	unclassified <- taxaAssg[unclId, ]  
	taxaAssg <- taxaAssg[-unclId, ]  
	
	# move the higher classification to the top of each taxa group and create taxagroup for colouring
	for ( taxaCol in (ncol(taxaAssg)-2):(colTotal+1) ) {
	    matchTaxaDF <- data.frame(row.names=unique(taxaAssg[,taxaCol]))
		matchTaxaDF$match1stcol <- match(rownames(matchTaxaDF),taxaAssg[,1])
		matchTaxaDF$matchtaxacol <- match(rownames(matchTaxaDF),taxaAssg[, taxaCol])
		# remove rows with NAs in data.frame
		matchTaxaDF <- matchTaxaDF[complete.cases(matchTaxaDF),]
				
		# category "high level" taxa higher than groupLevel
		if (nrow(matchTaxaDF) > 0) {
			if (taxaCol > colTaxaGroupLevel) {
				taxaAssg[matchTaxaDF[,1],"taxagroup"] <- highLevel
	#			taxaAssg[matchTaxaDF[,1],"space"] <- ""		
			} else {
				# add characters to 1st column taxa to build the rank hierarchy in axis label 
				notMatched <- setdiff(1:nrow(taxaAssg), matchTaxaDF[,1]) 
				taxaAssg[notMatched,"space"] <- paste(taxaAssg[notMatched,"space"], "--", sep = "")		
			}
				
			# remove rows with same taxa match index for simplification
			matchTaxaDF <- matchTaxaDF[matchTaxaDF[,1]!=matchTaxaDF[,2],]
		
			if (nrow(matchTaxaDF) > 0) {
				for ( mt in 1:nrow(matchTaxaDF) ) {
					taxaAssg <- mvRow(taxaAssg, matchTaxaDF[mt,1], matchTaxaDF[mt,2])
				}
			}
		}
	}
	
	taxaAssg <- rbind(unclassified, taxaAssg)
		
	rownames(taxaAssg) <- 1:nrow(taxaAssg)

	#####  make "Others" category #####
	# make "Others" by percThr: percentage threshold of total reads in taxaAssg$Total
	if (percThr > 0) {
	    totalThr <- sum(taxaAssg$Total) * percThr
		Others <- colSums(taxaAssg[which(taxaAssg$Total<=totalThr),2:colTotal])
		Others <- c("Others",Others,rep("Others", ncol(taxaAssg)-colTotal-1),"")	
		taxaAssgId <- which(taxaAssg$Total>totalThr)
		# avoid error: invalid factor level, NA generated
		taxaAssg <- data.frame(lapply(taxaAssg, as.character), stringsAsFactors=FALSE)
		taxaAssg <- rbind(taxaAssg[taxaAssgId,], Others)		
	}

	# make "Others" by maxNTaxa
#	if (nrow(taxaAssg) > maxNTaxa) {
#		Others <- colSums(taxaAssg[-(1:maxNTaxa),-1])
#		Others <- c("Others",Others)
#		taxaAssg <- data.frame(lapply(taxaAssg, as.character), stringsAsFactors=FALSE)
#		taxaAssg <- rbind(taxaAssg[1:maxNTaxa,], Others)		
#	}

	#####  prepare data frame for chart #####
	print( paste( "Select", length(unique(taxaPaths[,rankLevel])), "taxa from the total of ", rankLevel, "+ high-level taxa", nrow(taxaAssg), 
			", selected taxa groups =", length(unique(taxaAssg$taxagroup)), 
			", total OTUs =", nrow(taxaAssgReads), ", total reads =",  sum(taxaAssgReads[,-ncol(taxaAssgReads)]) ) )  

	xlab <- xlab( paste(nrow(taxaAssg), " (of ", length(unique(taxaPaths[,rankLevel])), ") ",  
					rankLevel, " + high-level taxa (reads > ", percThr*100, "% of total) ", sep = "") )

	taxaAssg[,1] <- paste(taxaAssg[,1], taxaAssg[,"space"], sep = "")
	
	y_string = "reads"
	# hard code for: 2nd col to Total
	# gather(taxaAssg[,1:colTotal], plot, reads, 2:colTotal)
	taxaAssgPerSample <- eval(parse(text = paste("gather(taxaAssg[,1:colTotal], plot", y_string, "2:colTotal)", sep=", ")))

	if ( nrow(taxaAssgPerSample) != (colTotal-1)*nrow(taxaAssg) ) 
		stop( paste("taxaAssgPerSample rows", nrow(taxaAssgPerSample), "should be", (colTotal-2)*nrow(taxaAssg), "!") )

	# assign taxa group to taxaAssgPerSample, including total column 
    taxaAssgPerSample$taxagroup <- rep(taxaAssg$taxagroup, colTotal-1)

    print(paste(length(unique(taxaAssg$taxagroup)), "taxa group from data :"))
    print(unique(taxaAssg$taxagroup))

	# choose the colours for taxagroup based on the union of all taxagroup found across all the markers
	if (expId == 1) {
		taxaGroupUnion <- c()
		platteeUsedId <- c()
		taxaUsedId <- c()
		platteeId <- 1:length(unique(taxaAssg$taxagroup))
		newPlatteeId <- platteeId			
	} else {
		# match previous union of taxa groups, if NA then the group not used previously
		platteeId <- match(unique(taxaAssg$taxagroup), taxaGroupUnion)	
		taxaUsedId <- platteeId[!is.na(platteeId)]
		platteeId[!is.na(platteeId)] <- platteeUsedId[taxaUsedId]

		if (is.na(all(platteeId))) {
			newPlatteeId <- 1:length(myPalette)
			newPlatteeId <- newPlatteeId[! newPlatteeId %in% platteeUsedId[taxaUsedId]]	
			# reuse the rest of group ids not being used by multi-makers		
			newPlatteeId <- newPlatteeId[1:length(platteeId[is.na(platteeId)])] 
			platteeId[is.na(platteeId)] <- newPlatteeId[1:length(platteeId[is.na(platteeId)])]
		} else {
			newPlatteeId <- c()
		}
	}

	taxaGroupUnion <- union(taxaGroupUnion, unique(taxaAssg$taxagroup))
	platteeUsedId <- c(platteeUsedId, newPlatteeId)
	
	if ( length(taxaGroupUnion) != length(platteeUsedId) ) 
		stop( paste("length of platteeUsedId ", length(platteeUsedId), "! = length of taxaGroupUnion", length(taxaGroupUnion), "!") )
 

	print(paste(length(taxaUsedId), "taxa group ids used previously :"))
	print(taxaUsedId)

	print(paste(length(taxaGroupUnion), "taxa group in current union set :"))
	print(taxaGroupUnion)

	# taxaGroupId only choose new id matched in taxaGroupUnion when the group is reused by >1 markers 
    print(paste(length(platteeId), "taxa group and their colour palette ids :"))
    print(platteeId)

	###### phyla OTUs bar chart #####
	taxaAssgPerSample[,y_string] <- as.numeric(taxaAssgPerSample[,y_string])
	taxaAssgPerSample[,rankLevel] = factor(taxaAssgPerSample[,rankLevel], rev(unique(taxaAssgPerSample[,rankLevel])))
	taxaAssgPerSample$taxagroup = factor(taxaAssgPerSample$taxagroup, unique(taxaAssgPerSample$taxagroup))

	breaks=c(0.1, 1, 10, 100, 1000)
	if (max(taxaAssgPerSample[,y_string]) < 500) {
	  breaks=c(0.1, 1, 10, 100)
	} else if (max(taxaAssgPerSample[,y_string]) < 50) {
	  breaks=c(0.1, 1, 10)
    }

	pdf(paste(workingPath, "figures/taxa-sequences-", matrixNames[expId], "-", rankLevel, ".pdf", sep = ""), width=18, height=8)	

	print( ggplot(taxaAssgPerSample, aes(x = order, y = y_string, fill = taxagroup)) + geom_bar(stat = "identity", position = "identity") + 
		 scale_y_continuous(trans = mylog_trans(base=10, from=-0.5), labels = function (x) floor(x), 
					breaks=breaks, expand = c(0, 0)) +
		 coord_flip() + theme_bw() + facet_grid( ~ plot) + ylab("Read abundance") + xlab + 
		 theme(legend.position=c(0.4, 1.09), legend.direction="horizontal", plot.margin=unit(c(2.0,0.3,0.2,0.8), "cm"), 
				 panel.margin = unit(0.8, "lines"), axis.title.y=element_text(vjust=2.8)) + 
		 guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
		 scale_fill_manual(name = paste(groupLevel, ":"),values = myPalette[platteeId]) )

	invisible(dev.off()) 


	###### by OTUs #####

	# create matrix based on OTU appearance in each plot
	#communityMatrixOTUs <- communityMatrix[communityMatrix>0]


	# create taxa table based on the number of OTUs
	#taxaAssgOTUs <-
}
