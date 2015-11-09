# Create a figure: the number of reads assigned to a pre-defined rank level and also high-level taxa, and coloured by another given rank level.
# Input: 1) taxonomy table from MEGAN, 1st column is OTUs, the last column is taxonomy, the middles are samples (plots/subplots);
#
# Author: Andrew Dopheide, Walter Xie
# Accessed on 16 Oct 2015

ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus")
# "#CCCCCC", "#999999"
myPalette <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#AE22EC", "#FAF629",
               "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#E69F00", "#FF6666", "#56B4E9", "#009E73", "#F0E442","#26ED26",    
               "#0072B2", "#D55E00", "#CC79A7", "#E6AB02", "#33ff99", "#fb9a99", "#e31a1c", "#fdbf6f", "#cab2d6", "#e41a1c", 
               "#377eb8", "#ff9933", "#4daf4a", "#984ea3", "#cc0033", "#a65628", "#9999FF", "#f781bf", "#a6cee3", "#1f78b4", 
               "#b2df8a", "#a6cee3", "#ffff33", "#006600", "#99AC06", "#2F3C56", "#A229FA", "#EC22DE", "#1DD3E8", "#9933FF")
total_string = "Total"

# rankLevel, groupLevel must be lower case to match ranks
# taxaAssg: 1st col is rankLevel, 2nd-? is CM, last 2nd is Total, last is groupLevel
# plotTotal: plot total section or not
# return a list
taxonomyAssignment <- function(taxaAssg, rankLevel="phylum", groupLevel="kingdom", percThr=0, 
                               legend_nrow=1, y_string="reads", plotTotal=TRUE, title=NULL) {

  cat("\nSet percentage threshold percThr =", percThr, ", and the value in bar represents \'", y_string, "\'.\n") 
  cat("rankLevel = ", rankLevel, ", groupLevel = ", groupLevel, ", legend_nrow = ", legend_nrow, 
      ", plotTotal = ", plotTotal, ", title = ", title, ".\n")
  
  ###### Preprocessing ######
  if (nrow(taxaAssg[taxaAssg[ ,total_string] < 1,]) > 0) {
    # rm rows rowSums==0 
    cat("Warning: removing", nrow(taxaAssg[taxaAssg[ ,total_string] < 1,]), "rows whose total is 0 :\n")
    
    taxaAssg <- taxaAssg[taxaAssg[ ,total_string] > 0,]
  }
  
  group = unique(taxaAssg[,groupLevel]) 
  cat("Input non-zero taxaAssg with", nrow(taxaAssg), rankLevel, ", assigned to", 
      length(group), groupLevel, ", total", y_string, "=",  sum(taxaAssg[ ,total_string]), ".\n")  
  
  xlab <- xlab( paste(nrow(taxaAssg), rankLevel, "(or higher-level taxon)") )
  ylab <- ylab(y_string)		
  
  colTotal <- which(colnames(taxaAssg)==total_string)
  
  #####  make "Others" category #####
  # make "Others" by percThr: percentage threshold of total reads in taxaAssg[ ,total_string]
  if (percThr > 0) {
    totalThr <- sum(taxaAssg[ ,total_string]) * percThr
    Others <- colSums(taxaAssg[which(taxaAssg[ ,total_string]<=totalThr),2:colTotal])
    Others <- c("Others",Others,rep("Others", ncol(taxaAssg)-colTotal),"")	
    taxaAssgId <- which(taxaAssg[ ,total_string]>totalThr)
    
    cat("Filter out", nrow(taxaAssg)-length(taxaAssgId), "taxa to Others category whose total <=", 
        percThr, "of the total of whole matrix.\n")  
    
    xlab <- xlab( paste(length(taxaAssgId), " of ", nrow(taxaAssg), rankLevel, " (", tolower(y_string), 
                        " > ", percThr*100, "% of total) ", sep = "") )		
    
    # avoid error: invalid factor level, NA generated
    taxaAssg <- data.frame(lapply(taxaAssg, as.character), check.names=FALSE, stringsAsFactors=FALSE)
    taxaAssg <- rbind(taxaAssg[taxaAssgId,], Others)
    
    group = unique(taxaAssg[,groupLevel]) # include Other
    # -1 to consider "Others" category in the last row
    cat("Select", nrow(taxaAssg)-1, rankLevel, "by percThr = ", percThr, ", assigned to", length(group), groupLevel, 
        ", total", y_string, "=",  sum(taxaAssg[ ,total_string]), ".\n")  
  }
  
  # make "Others" by maxNTaxa
  #	if (nrow(taxaAssg) > maxNTaxa) {
  #		Others <- colSums(taxaAssg[-(1:maxNTaxa),-1])
  #		Others <- c("Others",Others)
  #		taxaAssg <- data.frame(lapply(taxaAssg, as.character), stringsAsFactors=FALSE)
  #		taxaAssg <- rbind(taxaAssg[1:maxNTaxa,], Others)		
  #	}
  
  #####  prepare data frame for chart #####
  lastCol = colTotal
  if (!plotTotal) 
    lastCol = colTotal - 1
  # hard code for: 2nd col to Total
  # gather(taxaAssg[,1:lastCol], plot, reads, 2:lastCol)
  taxaAssgPerSample <- eval(parse(text = paste("gather(taxaAssg[,1:lastCol], plot", y_string, "2:lastCol)", sep=", ")))
  
  if ( nrow(taxaAssgPerSample) != (lastCol-1)*nrow(taxaAssg) ) 
    stop( paste("taxaAssgPerSample rows", nrow(taxaAssgPerSample), "should be", (lastCol-2)*nrow(taxaAssg), "!") )
  
  # assign taxa group to taxaAssgPerSample, including total column 
  taxaAssgPerSample[,groupLevel] <- rep(taxaAssg[,groupLevel], lastCol-1)
  
  cat(length(group), "taxa group from data :\n")
  cat(paste(group, collapse = ", "),"\n")
  
  ###### bar chart #####
  taxaAssgPerSample[,y_string] <- as.numeric(taxaAssgPerSample[,y_string])
  taxaAssgPerSample[,rankLevel] = factor(taxaAssgPerSample[,rankLevel], rev(unique(taxaAssgPerSample[,rankLevel])))
  taxaAssgPerSample[,groupLevel] = factor(taxaAssgPerSample[,groupLevel], unique(taxaAssgPerSample[,groupLevel]))
  
  maxLabelLen=max(nchar( as.character(taxaAssg[,rankLevel]) )) 
  # log scale: 1 cannot be plot
  fill_string = groupLevel
  #	legend_nrow = ceiling(length(unique(taxaAssgPerSample[,groupLevel])) / 11) # max 10 legend per row
  
  if (max(taxaAssgPerSample[,y_string]) < 50) {
    breaks=c(0.1, 1, 10)
  } else if (max(taxaAssgPerSample[,y_string]) < 500) {
    breaks=c(0.1, 1, 10, 100)
  } else if (max(taxaAssgPerSample[,y_string]) < 5000) {
    breaks=c(0.1, 1, 10, 100, 1000)
  } else if (max(taxaAssgPerSample[,y_string]) < 50000) {
    breaks=c(0.1, 1, 10, 100, 1000, 10000)
  } else if (max(taxaAssgPerSample[,y_string]) < 500000) {
    breaks=c(0.1, 1, 10, 100, 1000, 10000, 100000)
  } else {
    breaks=c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000)
  }
  
  if (length(group) > length(myPalette)) 
    myPalette <- rainbow_hcl(length(group))
    
  p <- ggplot(taxaAssgPerSample, aes_string(x = rankLevel, y = y_string, fill = fill_string)) + 
    geom_bar(stat = "identity", position = "identity") + 
    scale_y_continuous(trans = mylog_trans(base=10, from=-0.3), labels = scientific_10, 
                       breaks=breaks, expand = c(0, 0)) +
    coord_flip() + theme_bw() + facet_grid( ~ plot) + ylab + xlab + 
    theme(legend.position="top", legend.direction="horizontal", plot.margin=unit(c(0.2,0.5,0.2,0.8), "cm"), 
          panel.margin = unit(0.8, "lines"), axis.title.y=element_text(vjust=2.8), axis.text.x = element_text(vjust=-0.1)) + 
    guides(fill=guide_legend(nrow=legend_nrow,byrow=TRUE)) +
    scale_fill_manual(name = paste(groupLevel, ":"),values = myPalette) 
  
  if (!is.null(title)) 		
    p <- p + ggtitle(paste(title))	   

  cat("Return nrow =", nrow(taxaAssg), ", ncol =", ncol(taxaAssg), ", ngroup =", length(group), 
      ", maxlablen =", maxLabelLen, ".\n") 
  # Return a list containing the filename
  list(
    nrow = nrow(taxaAssg),
    ncol = ncol(taxaAssg), # include 'taxa' and 'group' col
    group = unique(taxaAssg[,groupLevel]), # include Other
    maxLabelLen = maxLabelLen, # max length of characters in label
    plot = p # ggplot
  )
}


# Fix x-axis number format
scientific_10 <- function(x) {
  text=gsub("1e\\+00", "1", scientific_format()(x))
  text=gsub("1e\\+01", "10", text)
  text=gsub("1e\\+", "10^", text)
  parse(text=text)
}

# taxaAssgReads = CM + rankLevel + groupLevel
# aggregate taxaAssgReads by rankLevel and y_string = OTU abundance or Read abundance
getTaxaAssg <- function(taxaAssgReads, rankLevel, groupLevel, y_string="reads") {
  taxaAssgReads[,groupLevel] <- gsub("root|cellular organisms|No hits|Not assigned", "unclassified", taxaAssgReads[,groupLevel])
  
  #y_string <- "OTU abundance" # alpha 0
  if (length(grep("OTU",y_string)) == 0) {
    # equivalent to d(x,lev=alpha,q=0)
    taxaAssg <- aggregate(as.formula(paste(". ~", rankLevel, "+", groupLevel)), data=taxaAssgReads[,-1], function(x) sum(x>0))
    #		} else if (grep("Effective OTU",y_string)) {
    #y_string <- "Effective OTU abundance" # alpha 1	
    #TODO: not working:	taxaAssg <- aggregate(as.formula(paste(".", rankLevel, sep=" ~ ")), data= tA.tmp, function(x) d(t(x),lev="gamma",q=1))	
  } else {
    #y_string <- "Read abundance" or "Reads"
    #grep("Read",y_string) 
    taxaAssg <- aggregate(as.formula(paste(". ~", rankLevel, "+", groupLevel)), data=taxaAssgReads[,-1], FUN=sum)	
  }
  
  # add Total column
  taxaAssg[,total_string] <- rowSums(taxaAssg[,-1:-2]) 
  
  # move groupLevel to the last col
  taxaAssg <- taxaAssg[,c(1,3:ncol(taxaAssg),2)]
  
  cat("Aggregate", nrow(taxaAssg), rankLevel, "belonging to", length(unique(taxaAssg[,groupLevel])), 
      groupLevel, "by counting", y_string, ".\n")

  return(taxaAssg)  
}

