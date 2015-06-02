# https://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/
library(ggplot2)
library(scales)
library(reshape2)
library(plyr)

# change config below
workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

divs <- c("gamma0", "beta1")
subTitles <- c(expression(paste("(a) "^"0", D[gamma])),expression(paste("(b) "^"1", D[beta])))

for (i in 1:length(divs)) {
	inputT <- paste(workingPath, "prob-4plots-table-", divs[i], ".txt", sep = "")
	#inputT <- paste(workingPath, "prob-4plots-table-beta1.txt", sep = "")

	rank_table <- read.table(inputT, sep="\t", header=T, row.names=1, check.names=FALSE)
	#rank_table[is.na(rank_table)] <- 0

	rank_table$Plots <- rownames(rank_table)
	rank.m <- melt(rank_table)

	colnames(rank.m)[grep("^value$", colnames(rank.m))] <- "Probability" # replace "value"
	rank.m$variable <- gsub("invertebrates", "inverts", rank.m$variable, ignore.case = T)
	rank.m$Plots <- factor(rank.m$Plots, levels=rev(unique(rank.m$Plots)))
	rank.m$variable <- factor(rank.m$variable, levels=unique(rank.m$variable))

#	pdf(paste(workingPath, "figures/prob-4plots-", divs[i], ".pdf", sep = ""), width=6, height=6)

	mdsp <- ggplot(rank.m, aes(x=variable, y=Plots, fill = Probability)) + 
			 geom_tile(colour = "white") + 
			 scale_fill_gradient(high = "steelblue", low = "white") + guides(fill = guide_legend(reverse=TRUE)) + ggtitle(subTitles[i]) +
			 theme(legend.position="top", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) 
	
	if (i==1) mylegend<-g_legend(mdsp) 		 
 
	mdsp <- mdsp + theme(legend.position="none")
	assign(paste('mdsp', i, sep=''), mdsp)
#	invisible(dev.off())         
}

pdf(paste(workingPath, "figures/prob-4plots-gamma0-beta1.pdf", sep = ""), width=8, height=5)	   
grid.arrange(mylegend,arrangeGrob(mdsp1,mdsp2,ncol = 2, nrow=1), ncol=1, nrow=2, heights=c(1/15,14/15))
invisible(dev.off()) 



