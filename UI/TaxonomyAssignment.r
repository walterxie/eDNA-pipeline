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

# http://stackoverflow.com/questions/22295253/force-bars-to-start-from-a-lower-value-than-0-in-ggplot-geom-bar-in-r
# defining the scale change
# scale_y_continuous(trans = mylog_trans(base=10, from=-2)) # starts from 1e-2
library(scales)
mylog_trans <- function(base=exp(1), from=0) {
  trans <- function(x) log(x, base)-from
  inv <- function(x) base^(x+from)
  trans_new("mylog", trans, inv, log_breaks(base=base), domain = c(base^from, Inf))
}


# how taxa to display, choose either maxNTaxa or percThr
#maxNTaxa=30 
#percThr = 0.001

# MEGAN mapping file column order
# no domain of "prokaryote", Bacteria's parent is "cellular organisms", and same as Archaea
#ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus")
rankLevel="Taxa" 
groupLevel="Group" # gives colour, and must higher than rankLevel
ylabTxt="Species abundance"

# "#CCCCCC", "#999999", 
myPalette <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", 
			   "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#E69F00", "#FF6666", "#56B4E9", "#009E73", "#F0E442", "#9933FF",   
				"#0072B2", "#D55E00", "#CC79A7", "#E6AB02", "#33ff99", "#fb9a99", "#e31a1c", "#fdbf6f", "#cab2d6", "#e41a1c", 
				"#377eb8", "#ff9933", "#4daf4a", "#984ea3", "#cc0033", "#a65628", "#9999FF", "#f781bf", "#a6cee3", "#1f78b4", 
				"#b2df8a", "#a6cee3", "#ffff33", "#006600")

#taxaGroupUnion <- c()
#, rankLevel=rankLevel, groupLevel=groupLevel, ylabTxt=ylabTxt
TaxonomyAssignment <- function(taxaAssg=taxaAssg, percThr=0.001, plotTitle=plotTitle, legend_nrow=1) {
  cat("Set percentage threshold percThr =", percThr, ".\n") 
  
  # rm white space
	colnames(taxaAssg) <- gsub(" ", "", colnames(taxaAssg), ignore.case = T)
	
	colnames(taxaAssg)[1] <- rankLevel
	colnames(taxaAssg)[ncol(taxaAssg)] <- groupLevel

	taxaAssg$Total <- rowSums(taxaAssg[,-c(1,ncol(taxaAssg))]) 
	
	# mv groupLevel col to the last col, easy to make "Others" category 
	colGroupLevel <- which(colnames(taxaAssg)==groupLevel)
	taxaAssg <- taxaAssg[,c(1:(colGroupLevel-1),ncol(taxaAssg),colGroupLevel)]
	
	# rm rows rowSums==0 
	cat("Remove rows whose total is 0 :\n")
	print(taxaAssg[taxaAssg$Total < 1,])
	
	taxaAssg <- taxaAssg[taxaAssg$Total > 0,]
	
	cat("Input non-zero community matrix with", nrow(taxaAssg), "taxa, assigned to", length(unique(taxaAssg[,groupLevel])), 
			"taxa groups, total individuals/reads =",  sum(taxaAssg$Total), ".\n")  
		
	xlab <- xlab( paste(nrow(taxaAssg), " taxa", sep = "") )
	ylab <- ylab(ylabTxt)		

	colTotal <- which(colnames(taxaAssg)=="Total")
	    
	#####  make "Others" category #####
	# make "Others" by percThr: percentage threshold of total reads in taxaAssg$Total
	if (percThr > 0) {
		totalThr <- sum(taxaAssg$Total) * percThr
		Others <- colSums(taxaAssg[which(taxaAssg$Total<=totalThr),2:colTotal])
		Others <- c("Others",Others,rep("Others", ncol(taxaAssg)-colTotal),"")	
		taxaAssgId <- which(taxaAssg$Total>totalThr)

		cat( "Filter out", nrow(taxaAssg)-length(taxaAssgId), "taxa to Others category whose total <=", percThr, "of the total of whole matrix.\n")  

		xlab <- xlab( paste(length(taxaAssgId), " of ", nrow(taxaAssg), " taxa (", tolower(ylabTxt), " > ", percThr*100, "% of total) ", sep = "") )		
		
		# avoid error: invalid factor level, NA generated
		taxaAssg <- data.frame(lapply(taxaAssg, as.character), check.names=FALSE, stringsAsFactors=FALSE)
		taxaAssg <- rbind(taxaAssg[taxaAssgId,], Others)
		
		# -1 to consider "Others" category in the last row
		cat("Select", nrow(taxaAssg)-1, "taxa, assigned to", length(unique(taxaAssg[,groupLevel]))-1, 
			"taxa groups, total individuals/reads =",  sum(as.numeric(taxaAssg$Total)), ".\n")  
	}

	# make "Others" by maxNTaxa
#	if (nrow(taxaAssg) > maxNTaxa) {
#		Others <- colSums(taxaAssg[-(1:maxNTaxa),-1])
#		Others <- c("Others",Others)
#		taxaAssg <- data.frame(lapply(taxaAssg, as.character), stringsAsFactors=FALSE)
#		taxaAssg <- rbind(taxaAssg[1:maxNTaxa,], Others)		
#	}

	#####  prepare data frame for chart #####
	y_string = "reads"
	# hard code for: 2nd col to Total
	# gather(taxaAssg[,1:colTotal], plot, reads, 2:colTotal)
	taxaAssgPerSample <- eval(parse(text = paste("gather(taxaAssg[,1:colTotal], plot", y_string, "2:colTotal)", sep=", ")))
	
	if ( nrow(taxaAssgPerSample) != (colTotal-1)*nrow(taxaAssg) ) 
		stop( paste("taxaAssgPerSample rows", nrow(taxaAssgPerSample), "should be", (colTotal-2)*nrow(taxaAssg), "!") )

	# assign taxa group to taxaAssgPerSample, including total column 
	taxaAssgPerSample[,groupLevel] <- rep(taxaAssg[,groupLevel], colTotal-1)

	cat(length(unique(taxaAssg[,groupLevel])), "taxa group from data :\n")
	cat(paste(unique(taxaAssg[,groupLevel]), collapse = ", "),"\n")

	###### phyla OTUs bar chart #####
	taxaAssgPerSample[,y_string] <- as.numeric(taxaAssgPerSample[,y_string])
	taxaAssgPerSample[,rankLevel] = factor(taxaAssgPerSample[,rankLevel], rev(unique(taxaAssgPerSample[,rankLevel])))
	taxaAssgPerSample[,groupLevel] = factor(taxaAssgPerSample[,groupLevel], unique(taxaAssgPerSample[,groupLevel]))

#	pdfHeight=2+nrow(taxaAssg)*0.12
# maxLabelLen=max(nchar(taxaAssg[,rankLevel])) 
#	pdfWidth=0.1+maxLabelLen/10+(ncol(taxaAssg)-2)*1.5 
	
#	pdf(paste("figures/", "taxa-", plotTitle, ".pdf", sep = ""), width=pdfWidth, height=pdfHeight)	

	# log scale: 1 cannot be plot
	fill_string = groupLevel
#	legend_nrow = ceiling(length(unique(taxaAssgPerSample[,groupLevel])) / 11) # max 10 legend per row

	breaks=c(0.1, 1, 10, 100, 1000)
	if (max(taxaAssgPerSample[,y_string]) < 500) {
	  breaks=c(0.1, 1, 10, 100)
	} else if (max(taxaAssgPerSample[,y_string]) < 50) {
	  breaks=c(0.1, 1, 10)
	}
	
	p <- ggplot(taxaAssgPerSample, aes_string(x = rankLevel, y = y_string, fill = fill_string)) + 
			geom_bar(stat = "identity", position = "identity") + 
	    scale_y_continuous(trans = mylog_trans(base=10, from=-0.3), labels = function (x) floor(x), 
	                     breaks=breaks, expand = c(0, 0)) +
	    coord_flip() + theme_bw() + facet_grid( ~ plot) + ylab + xlab + 
#			ggtitle(paste(plotTitle)) +
			theme(legend.position="top", legend.direction="horizontal", plot.margin=unit(c(0.2,0.5,0.2,0.8), "cm"), 
				panel.margin = unit(0.8, "lines"), axis.title.y=element_text(vjust=2.8)) + 
			guides(fill=guide_legend(nrow=legend_nrow,byrow=TRUE)) +
			scale_fill_manual(name = paste(groupLevel, ":"),values = myPalette) 
#	print(p)
#	invisible(dev.off()) 
}
