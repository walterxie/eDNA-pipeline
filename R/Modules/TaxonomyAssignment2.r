# Create a figure: the number of reads assigned to a pre-defined rank level and also high-level taxa, 
# and coloured by another given rank level.
#
# Author: Andrew Dopheide, Walter Xie
# Accessed on 10 Nov 2015

ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus")
# "#CCCCCC", "#999999"
myPalette <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#AE22EC", "#FAF629",
               "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#E69F00", "#FF6666", "#56B4E9", "#009E73", "#F0E442","#26ED26",    
               "#0072B2", "#D55E00", "#CC79A7", "#E6AB02", "#33ff99", "#fb9a99", "#e31a1c", "#fdbf6f", "#cab2d6", "#e41a1c", 
               "#377eb8", "#ff9933", "#4daf4a", "#984ea3", "#cc0033", "#a65628", "#9999FF", "#f781bf", "#a6cee3", "#1f78b4", 
               "#b2df8a", "#a6cee3", "#ffff33", "#006600", "#99AC06", "#2F3C56", "#A229FA", "#EC22DE", "#1DD3E8", "#9933FF")
total_string = "Total"
# named by melt function
y_string="value"

tax_ref <- getTaxaRef()

# plot OTUs and reads in one graph, use column "var" to store flag
# rankLevel: the rank chosen to melt down community matrix
# groupLevel: gives colour, and must higher than rankLevel
# rankLevel, groupLevel must be lower case to match ranks
# taxaAssg: 1st col is rankLevel, 2nd-? is CM, last 3nd is Total, last 2nd is groupLevel, last is var (= OTUs/reads)
# plotTotal: plot total section or not
taxonomyAssignment2 <- function(taxaAssg, rankLevel="phylum", groupLevel="kingdom", percThr=0, 
                               legend_nrow=1, plotTotal=TRUE, title=NULL) {

  cat("\nSet percentage threshold percThr =", percThr, ", rankLevel = ", rankLevel, ", groupLevel = ", groupLevel, 
      ", legend_nrow = ", legend_nrow, ", plotTotal = ", plotTotal, ", title = ", title, ".\n")
  
  ###### Preprocessing ######
  if (nrow(taxaAssg[taxaAssg[ ,total_string] < 1,]) > 0) {
    # rm rows rowSums==0 
    cat("Warning: removing", nrow(taxaAssg[taxaAssg[ ,total_string] < 1,]), "rows whose total is 0 :\n")
    
    taxaAssg <- taxaAssg[taxaAssg[ ,total_string] > 0,]
  }
  
  taxa = unique(taxaAssg[,rankLevel]) 
  group = unique(taxaAssg[,groupLevel]) 
  cat("Input non-zero taxaAssg with", length(taxa), rankLevel, ", assigned to", length(group), groupLevel, ".\n")  
  
  xlab <- xlab( paste(length(taxa), rankLevel, "(or higher-level taxon)") )
  ylab <- ylab(paste("OTUs or sequences"))		

  #####  prepare data frame for chart #####
  if (!plotTotal) 
    taxaAssg <- taxaAssg[,-which(names(taxaAssg) %in% total_string)]
  
  taxaAssgPerSample <- melt(taxaAssg, id=c(rankLevel,groupLevel,"var"))
  taxaAssgPerSample[,rankLevel] <- factor(taxaAssgPerSample[,rankLevel], ordered = TRUE, 
                                          levels = rev(unique(tax_ref[,rankLevel])))
  taxaAssgPerSample[,groupLevel] <- factor(taxaAssgPerSample[,groupLevel], ordered = TRUE, 
                                           levels = unique(tax_ref[,groupLevel]))
  
  ###### bar chart #####
  maxLabelLen=max(nchar( as.character(taxaAssg[,rankLevel]) )) 

  if (length(group) > length(myPalette)) 
    myPalette <- rainbow_hcl(length(group))
  
  inter <- paste0("interaction(", rankLevel, ", variable)")
  
  p <- ggplot(taxaAssgPerSample) + 
    geom_point(data = taxaAssgPerSample[taxaAssgPerSample$var == "reads",], 
               aes_string(x = rankLevel, y = y_string, colour = groupLevel), shape = 5, size = 2) +
    geom_point(data = taxaAssgPerSample[taxaAssgPerSample$var == "OTUs",], 
               aes_string(x = rankLevel, y = y_string, colour = groupLevel), shape = 1, size = 2) +
    #geom_point(data = taxaAssgPerSample[var == "reads_min1"], aes(x = Phylum, y = value, colour = Kingdom), shape = 17, size = 2) +
    #geom_point(data = taxaAssgPerSample[var == "OTUs_min1"], aes(x = Phylum, y = value, colour = Kingdom), shape = 16, size = 2) +
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), label=scientific_10) + 
    facet_grid( ~ variable) + coord_flip()+ theme_bw() + xlab + ylab +
    theme(legend.position="top", legend.direction="horizontal",
          plot.margin=unit(c(0.2,0.5,0.2,0.8), "cm"), panel.margin = unit(0.8, "lines"), 
          axis.title.y=element_text(vjust=2), axis.title.x=element_text(vjust=-0.2), 
          axis.text.x = element_text(vjust=-0.1)) +
    geom_line(data = taxaAssgPerSample, aes_string(x = rankLevel, y = y_string, colour = groupLevel,
                                                   group=inter), size = 0.5, alpha = 0.5) 
  
  if (!is.null(title)) 		
    p <- p + ggtitle(paste(title))	   

  cat("Return nrow =", length(taxa), ", ncol =", ncol(taxaAssg), ", ngroup =", length(group), 
      ", maxlablen =", maxLabelLen, ".\n") 
  # Return a list containing the filename
  list(
    nrow = length(taxa), # number of taxa
    ncol = ncol(taxaAssg), # include 'taxa' and 'group' col
    group = unique(taxaAssg[,groupLevel]), # include Other
    maxLabelLen = maxLabelLen, # max length of characters in label
    plot = p # ggplot
  )
}

# taxaAssgReads = CM + rankLevel + groupLevel
# aggregate taxaAssgReads by rankLevel 
getTaxaAssg2 <- function(taxaAssgReads, rankLevel, groupLevel) {
  taxaAssgReads[,groupLevel] <- gsub("root|cellular organisms|No hits|Not assigned", "unclassified", taxaAssgReads[,groupLevel])
  
  # remove 1st column "Row.names" and last colmun groupLevel
  otus <- aggregate(as.formula(paste(". ~", rankLevel, "+", groupLevel)), data=taxaAssgReads[,-1], function(x) sum(x>0))
  reads <- aggregate(as.formula(paste(". ~", rankLevel, "+", groupLevel)), data=taxaAssgReads[,-1], FUN=sum)	
  
  # add Total column
  otus[,total_string] <- rowSums(otus[,-1:-2]) 
  reads[,total_string] <- rowSums(reads[,-1:-2]) 
  
  # move groupLevel to the last col
  otus <- otus[,c(1,3:ncol(otus),2)]
  reads <- reads[,c(1,3:ncol(otus),2)]
  otus$var <- "OTUs"
  reads$var <- "reads"
  
  cat("Aggregate", nrow(otus), rankLevel, "belonging to", length(unique(otus[,groupLevel])), groupLevel, ".\n")
  
  return(rbind(otus,reads))  
}

