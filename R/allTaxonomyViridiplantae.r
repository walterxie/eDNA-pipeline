# Author: Walter Xie
# Accessed on 23 June 2015

library(gplots)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(tidyr)
library(vegetarian)

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

# insert row to data frame
insertRow <- function(rowDF, to, taxaAssg) {	
	if (to == 1) {
		return(rbind(rowDF,taxaAssg))
	} else if (to > nrow(taxaAssg)) {
		return(rbind(taxaAssg,rowDF))
	} else {
		return( rbind(rbind(taxaAssg[1:(to-1),],rowDF), taxaAssg[to:(nrow(taxaAssg)),]) )
	}
}

# http://stackoverflow.com/questions/22295253/force-bars-to-start-from-a-lower-value-than-0-in-ggplot-geom-bar-in-r
# defining the scale change
# scale_y_continuous(trans = mylog_trans(base=10, from=-2)) # starts from 1e-2
library(scales)
mylog_trans <- function(base=exp(1), from=0) {
  trans <- function(x) log(x, base)-from
  inv <- function(x) base^(x+from)
  trans_new("mylog", trans, inv, log_breaks(base=base), domain = c(base^from, Inf))
}

# change config below
workingPath <- "~/svn/compevol/research/NZGenomicObservatory/MiSeq/DOCt1/"
matrixNames <-  c("18S","26S") # only for cm file name and folder name   
taxaFiles <- c("DOCx9_18S_PEARcontigs_ME1.0_OTUs_num_nt_taxatable.txt","DOCx9_26S_R1_reads_OTUs_num_nt_taxatable.txt")
taxaManualFiles <- c("LUCAS_18S_Viridaeplantae.txt","LUCAS_26S_Viridaeplantae.txt")

# how taxa to display, choose either maxNTaxa or percThr
#maxNTaxa=30 
percThr = 0
# not make a plot for the data set whose reads <= minReads
minReads = 10

belongTo = "Viridiplantae" 

# MEGAN mapping file column order
# c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus")
ranks <- c("division", "order", "family", "genus","species")
highLevel="high level"
rankLevel="species" 
groupLevel="division" # gives colour, and must higher than rankLevel

displayEmptyParent=FALSE #TRUE

myPalette <- c("#999999", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", 
			   "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#E69F00", "#FF6666", "#56B4E9", "#009E73", "#F0E442", "#9933FF",   
				"#0072B2", "#D55E00", "#CC79A7", "#E6AB02", "#33ff99", "#fb9a99", "#e31a1c", "#fdbf6f", "#cab2d6", "#e41a1c", 
				"#377eb8", "#ff9933", "#4daf4a", "#984ea3", "#cc0033", "#a65628", "#9999FF", "#f781bf", "#a6cee3", "#1f78b4", 
				"#b2df8a", "#a6cee3", "#ffff33", "#006600")

taxaGroupUnion <- c()

ylabs <- c("Read abundance", "OTU abundance")#, "Effective OTU abundance")
filePrefixes <- c("read", "gamma0")#, "gamma1")

for (ylabId in 1:length(ylabs)) {
ylab <- ylab(ylabs[ylabId])
filePrefix = filePrefixes[ylabId]
print(paste("Analyse", ylab, ", and give file prefix : ", filePrefix, "-"))

#for (belongId in 1:length(belongToList)) {
#	belongTo = belongToList[belongId]
#	rankLevel = rankLevelList[belongId]
#	groupLevel = groupLevelList[belongId] 
	cat(paste("\nChoose rank level", rankLevel, "grouped by", groupLevel, "for", belongTo, " :\n"))

	for (expId in 1:length(matrixNames)) {
		cat(paste("\nLoad data for", matrixNames[expId], " :\n"))
	
		##### load data #####
		inputCM <- paste(workingPath, "data/", matrixNames[expId], ".txt", sep="") 
		communityMatrix <- read.table(inputCM, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
		communityMatrix <- communityMatrix[order(rownames(communityMatrix)),]
		communityMatrix <- communityMatrix[,order(colnames(communityMatrix))]
		# U115 = U155 (correct name) and OJ95 = CJ95 (correct name). AL132 is the "difficult" site extracting DNA 
		colnames(communityMatrix) <- gsub("OJ95", "CJ95", colnames(communityMatrix), ignore.case = T)
		colnames(communityMatrix) <- gsub("U115", "U155", colnames(communityMatrix), ignore.case = T)
	
		# taxaPathsMEGAN only for validate rownames here
		inputTaxa <- paste(workingPath, "Taxonomy_tables/", taxaFiles[expId], sep="")
		taxaPathsMEGAN <- read.table(inputTaxa, header=TRUE, row.names=1, sep = "\t", stringsAsFactors=FALSE)  
		taxaPathsMEGAN <- taxaPathsMEGAN[order(rownames(taxaPathsMEGAN)),]

		rownames(communityMatrix) <- gsub(";size.*;", "", rownames(communityMatrix))
		rownames(taxaPathsMEGAN) <- gsub(";size.*;", "", rownames(taxaPathsMEGAN))

		##### filter out rows not belong to given taxa belongTo ##### 
		taxaPathsMEGAN <- taxaPathsMEGAN[which(grepl(belongTo, taxaPathsMEGAN[,1])),] # taxaPaths[,1] is taxa path separated by ;
	
		# taxaPaths here is manually classified
		inputTaxaManual <- paste(workingPath, "data/", taxaManualFiles[expId], sep="")
		taxaPaths <- read.table(inputTaxaManual, header=TRUE, row.names=1, sep = "\t", stringsAsFactors=FALSE)
		rownames(taxaPaths) <- paste("OTU_", rownames(taxaPaths), sep="") 
		taxaPaths <- taxaPaths[order(rownames(taxaPaths)),]
		colnames(taxaPaths) <- tolower(colnames(taxaPaths))
		
		# taxaPaths is a subset of taxaPathsMEGAN
		if (nrow(taxaPaths) == nrow(taxaPathsMEGAN)) {
			if ( all( tolower(rownames(taxaPaths)) != tolower(rownames(taxaPathsMEGAN)) ) ) 
				stop( paste("OTU names in", matrixNames[expId], "two taxa path files are not matched !") )
		} else if (nrow(taxaPaths) > nrow(taxaPathsMEGAN)) {
			stop( paste("Manual classification file cannot have more OTUs than MEGAN file in", matrixNames[expId], " !") )
		} else {
			if ( any( is.na( match(rownames(taxaPaths),rownames(taxaPathsMEGAN)) ) ) ) 
				stop( paste("OTU names in", matrixNames[expId], "two taxa path files are not matched !") )		
		}	

		# replace empty values with value from the higher rank column
		for (cl in 2:ncol(taxaPaths)) {
			taxaPaths[ taxaPaths[,cl] == "", cl  ]  <- taxaPaths[ taxaPaths[,cl] == "", cl-1 ]
		}

		# classified by Viridiplantae from Genbank, but corrected as Fungi by Ben
		numFungi <- nrow(taxaPaths[ taxaPaths[,1] == "Fungi", ])
		print( paste(numFungi, "OTUs are classified by Viridiplantae from Genbank, but corrected as Fungi by Ben.") )
		
		# remove all Fungi
		taxaPaths <- taxaPaths[ taxaPaths[,1] != "Fungi", ]
	
		# not make a plot for the data set whose reads < minReads
		if (nrow(taxaPaths) < minReads) {
			print( paste("Taxa of", matrixNames[expId], "NOT have classification for", belongTo, ", jump to next in the loop.") )	
			next
		}

		communityMatrix <- communityMatrix[match(rownames(taxaPaths),rownames(communityMatrix)),]

		if ( all( tolower(rownames(communityMatrix)) != tolower(rownames(taxaPaths)) ) ) 
			stop( paste("OTU names in", matrixNames[expId], "community matrix file and taxa path file are not matched !") )

		if ( ! rankLevel %in% colnames(taxaPaths) ) 
			stop( paste("Column name", rankLevel, "not exist in taxa path file for ", matrixNames[expId]) )
		if (! groupLevel %in% colnames(taxaPaths) ) 
			stop( paste("Column name", groupLevel, "not exist in taxa path file for ", matrixNames[expId]) )
	
		###### taxa assignment by reads #####
		# create taxa table based on the number of reads
		taxaAssgReads <- cbind(communityMatrix, taxaPaths[,rankLevel])
		colnames(taxaAssgReads)[ncol(taxaAssgReads)] <- rankLevel

		if (ylab=="OTU abundance") {
			#ylab <- "OTU abundance" # alpha 0
			# equivalent to d(x,lev=alpha,q=0)
			taxaAssg <- aggregate(as.formula(paste(".", rankLevel, sep=" ~ ")), data= taxaAssgReads, function(x) sum(x>0))
#		} else if (ylab=="Effective OTU abundance") {
			#ylab <- "Effective OTU abundance" # alpha 1	
#TODO: not working:	taxaAssg <- aggregate(as.formula(paste(".", rankLevel, sep=" ~ ")), data= taxaAssgReads, function(x) d(t(x),lev="gamma",q=1))	
		} else {
			#ylab <- "Read abundance" 
			taxaAssg <- aggregate(as.formula(paste(".", rankLevel, sep=" ~ ")), data= taxaAssgReads, FUN=sum)		
		}
		taxaAssg$Total <- rowSums(taxaAssg[,-1]) 

		# colnames(taxaPaths): "division" "order"    "family"   "genus"    "species"
		colRanks <- match(ranks, colnames(taxaPaths))
		colRanks <- colRanks[!is.na(colRanks)] 
		colRankLevel <- which(colnames(taxaPaths)==rankLevel)
		colGroupLevel <- which(colnames(taxaPaths)==groupLevel)
	
		if (colRankLevel <= colGroupLevel || colGroupLevel < colRanks[1]) 
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
	
		#### move the higher classification to the top of each taxa group and create taxagroup for colouring
		for ( taxaCol in (ncol(taxaAssg)-2):(colTotal+1) ) {
			matchTaxaDF <- data.frame(row.names=unique(taxaAssg[,taxaCol]))
			matchTaxaDF$match1stcol <- match(rownames(matchTaxaDF),taxaAssg[,1])
			matchTaxaDF$matchtaxacol <- match(rownames(matchTaxaDF),taxaAssg[, taxaCol])
			
			# get empty parent
			emptyParentDF <- matchTaxaDF[is.na(matchTaxaDF[,1]),]
			# remove rows with NAs in data.frame
			matchTaxaDF <- matchTaxaDF[complete.cases(matchTaxaDF),]
				
			# category "high level" taxa higher than groupLevel
			if (nrow(matchTaxaDF) > 0) {
				if (taxaCol > colTaxaGroupLevel) {
					taxaAssg[matchTaxaDF[,1],"taxagroup"] <- highLevel
				} else {
					# add characters to 1st column taxa to build the rank hierarchy in axis label 
					notMatched <- setdiff(1:nrow(taxaAssg), matchTaxaDF[,1]) 
				}
				
				# remove rows with same taxa match index for simplification
				matchTaxaDF <- matchTaxaDF[matchTaxaDF[,1]!=matchTaxaDF[,2],]
		
				if (nrow(matchTaxaDF) > 0) {
					for ( mt in 1:nrow(matchTaxaDF) ) {
						taxaAssg <- mvRow(taxaAssg, matchTaxaDF[mt,1], matchTaxaDF[mt,2])
					}
				}
			}
			
			# insert empty parent
			if (displayEmptyParent && nrow(emptyParentDF) > 0) {
				for ( insRow in 1:nrow(emptyParentDF) ) {
					to <- emptyParentDF[insRow,2]+insRow-1
					rowDF <- taxaAssg[to,]
					# hard code
					rowDF[2:colTotal]=0
					rowDF[1]=rowDF[taxaCol]
					rowDF[(colTotal+1):taxaCol]=rowDF[taxaCol]
					taxaAssg <- insertRow(rowDF, to, taxaAssg)			
				}			
			}
			
		}
		
		rownames(taxaAssg) <- 1:nrow(taxaAssg)
		
		#### add space to "space" column
		space <- " -"
		for ( taxaCol in colTaxaGroupLevel:(colTotal+1) ) {
			if (taxaCol == (colTotal+1)) {
				notMatched <- which(taxaAssg[,1]!=taxaAssg[,colTotal+1]) 
			} else {
				notMatched <- which(taxaAssg[,taxaCol]!=taxaAssg[,taxaCol-1]) 
			}
			taxaAssg[notMatched,"space"] <- paste(taxaAssg[notMatched,"space"], space, sep = "")	
		}

		##### add genus to species name regarding taxonomy rule #####
		if ( all( is.element(c("genus","species"), tolower(colnames(taxaAssg))) ) ) {
			colSpecies <- which(colnames(taxaAssg)=="species")
			colGenus <- which(colnames(taxaAssg)=="genus")
			# factor causing error 
			taxaAssg[, colSpecies] <- sapply(taxaAssg[, colSpecies], as.character)
			
			isSpecies <- which(taxaAssg[,colSpecies] != taxaAssg[,colGenus])
			# keep at most 10 characters from genus
			prefix <- substring(taxaAssg[isSpecies, colGenus], 1, 10)
			taxaAssg[isSpecies, colSpecies] <- paste(prefix, taxaAssg[isSpecies, colSpecies], sep=".")
		}
	

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
		print( paste( "Select", length(unique(taxaPaths[,rankLevel])), "taxa from the total of", rankLevel, "+ high-level taxa", nrow(taxaAssg), 
				", selected taxa groups =", length(unique(taxaAssg$taxagroup)), 
				", total OTUs =", nrow(taxaAssgReads), ", total reads =",  sum(taxaAssgReads[,-ncol(taxaAssgReads)]) ) )  

		if (percThr > 0) {
			xlab <- xlab( paste(nrow(taxaAssg), " (of ", length(unique(taxaPaths[,rankLevel])), ") ",  
							rankLevel, " + high-level taxa (reads > ", percThr*100, "% of total) ", sep = "") )
		} else {
			if (belongTo == "Metazoa") { # non Arthropoda
				xlab <- xlab( paste(nrow(taxaAssg), " ", rankLevel, " + high-level taxa (excluding Arthropoda)", sep = "") )
			} else {
				xlab <- xlab( paste(nrow(taxaAssg), " ", rankLevel, " + high-level taxa", sep = "") )
			}
		}

		taxaAssg[,1] <- paste(taxaAssg[,1], taxaAssg[,"space"], sep = "")

		# MEAGEN bug: show Eukaryota in all rank columns for Viridiplantae, whose taxa path contains Viridiplantae 
		taxaAssg[,1] <- gsub("Eukaryota", "Viridiplantae", taxaAssg[,1])

		# hard code for: 2nd col to 11th, or AU127:Total
		taxaAssgPerSample <- gather(taxaAssg[,1:colTotal], plot, reads, 2:colTotal)

		if ( nrow(taxaAssgPerSample) != (colTotal-1)*nrow(taxaAssg) ) 
			stop( paste("taxaAssgPerSample rows", nrow(taxaAssgPerSample), "should be", (colTotal-2)*nrow(taxaAssg), "!") )

		# assign taxa group to taxaAssgPerSample, including total column 
		taxaAssgPerSample$taxagroup <- rep(taxaAssg$taxagroup, colTotal-1)

		print(paste(length(unique(taxaAssg$taxagroup)), "taxa group from data :"))
		print(unique(taxaAssg$taxagroup))

		# remove grey if no highLevel
		if (!is.element(highLevel, taxaAssg[,"taxagroup"])) 
			myPalette <- myPalette [! myPalette %in% "#999999"]

		# choose the colours for taxagroup based on the union of all taxagroup found across all the markers
		if (length(taxaGroupUnion)==0) {
			platteeUsedId <- c()
			taxaUsedId <- c()
			# also use grey for "high level"
	#		hLId <- which(unique(taxaAssg$taxagroup)==highLevel)
	#		if (hLId == 1) {
				platteeId <- 1:length(unique(taxaAssg$taxagroup))
	#		} else if (hLId == length(unique(taxaAssg$taxagroup))) {
	#			platteeId <- c(2:hLId,1)
	#		} else {
	#			platteeId <- c(2:hLId,1,(hLId+1):length(unique(taxaAssg$taxagroup)))	
	#		}		
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
		taxaAssgPerSample$reads <- as.numeric(taxaAssgPerSample$reads)
		taxaAssgPerSample[,rankLevel] = factor(taxaAssgPerSample[,rankLevel], rev(unique(taxaAssgPerSample[,rankLevel])))
		taxaAssgPerSample$taxagroup = factor(taxaAssgPerSample$taxagroup, unique(taxaAssgPerSample$taxagroup))

		pdfHeight=2+nrow(taxaAssg)*0.12
		pdf(paste(workingPath, "figures/", filePrefix, "-", belongTo, "-", matrixNames[expId], "-", rankLevel, "-manual.pdf", sep = ""), 
			width=16, height=pdfHeight)	

		# log scale: 1 cannot be plot
		y_string = "reads"
		fill_string = "taxagroup"
		legend_nrow = ceiling(length(unique(taxaAssgPerSample$taxagroup)) / 11) # max 10 legend per row
		if (belongTo=="Bacteria") 
			legend_nrow = ceiling(length(unique(taxaAssgPerSample$taxagroup)) / 7) # max 6 legend per row

		print( ggplot(taxaAssgPerSample, aes_string(x = rankLevel, y = y_string, fill = fill_string)) + 
				geom_bar(stat = "identity", position = "identity") + 
				scale_y_continuous(trans = mylog_trans(base=10, from=-0.25), labels = function (x) floor(x), 
					breaks=c(0.1, 1, 10, 100), expand = c(0, 0)) +
				coord_flip() + theme_bw() + facet_grid( ~ plot) + ylab + xlab + 
				ggtitle(paste(matrixNames[expId], "data set for", belongTo)) +
				theme(legend.position="top", legend.direction="horizontal", plot.margin=unit(c(0.2,0.5,0.2,0.8), "cm"), 
					panel.margin = unit(0.8, "lines"), axis.title.y=element_text(vjust=2.8)) + 
				guides(fill=guide_legend(nrow=legend_nrow,byrow=TRUE)) +
				scale_fill_manual(name = paste(groupLevel, ":"),values = myPalette[platteeId]) )

		invisible(dev.off()) 

	} #END for expId
#} #END for belongId
} #END for ylabelId

