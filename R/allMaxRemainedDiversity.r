
library(gplots)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(xtable)


if(!exists("tableFile")) stop("table file is missing !")
if(!exists("figDir")) stop("figure folder name is missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("otuThr")) otuThr = 97

#source("Modules/init.R", local=TRUE)
source("Modules/MaximizeDiversity.R")

env <- getSampleMetaData(TRUE)
env[,"ForestType"] <- gsub(":.*", "", env[,"ForestType"], ignore.case = T)
env[,"ForestType"] <- gsub("x", "unknown", env[,"ForestType"], ignore.case = T)

######## heatmap #######
divs <- c("gamma-0", "beta-1")
subTitles <- c(expression(paste("Ranked by "^"0", D[gamma])),expression(paste("Ranked by "^"1", D[beta])))

taxa.group="assigned"
for (i in 1:length(divs)) {
  # ranks, last 2 col are mean and sd
  allRanks <- getMaxRemainedDiversityRank(divs[i], taxa.group) 
  rank.elv <- merge(allRanks, env[,c("Elevation","ForestType")], by = "row.names")
  rank.elv <- rank.elv[order(rank.elv[,"Elevation"]),]
  colnames(rank.elv)[1] <- "Samples"

  fname <- paste("max-remained", divs[i], taxa.group,"heatmap.pdf", sep = "-")
  plotHeatmapMaxRemDivRank(rank.elv, fname, title=subTitles[i])
  
  rownames(rank.elv) <- rank.elv[,"Samples"]
  rank.elv <- rank.elv[,-c(1,ncol(rank.elv)-1,ncol(rank.elv))]
  printXTable(rank.elv, caption = paste("Ranking sampling plots by removing plots sequentially so as to 
			 minimize the loss of overall", divs[i], "diversity among the remaining plots. 
	     1 is the most important and removed at the last, 10 is the least important and removed in the beginning."), 
       label = paste0("tab:maxRemained", divs[i]), file=tableFile)
}

taxaGroups <- c("ANIMALIA", "FUNGI", "PROTISTS")
for (i in 1:length(divs)) {
  # hard code
  allRanks <- getMaxRemainedDiversityRank(divs[i], taxa.group="BACTERIA") 
  allRanksGourp <- data.frame(row.names = rownames(allRanks), stringsAsFactors = F)
  allRanksGourp[,"16S.BACTERIA"] <- allRanks[,"16S"]
  
  for (taxag in taxaGroups) {
      taxa.group <- taxag
      allRanks <- getMaxRemainedDiversityRank(divs[i], taxa.group=taxa.group) 
      colnames(allRanks) <- paste(colnames(allRanks), taxa.group, sep = ".")
      
      allRanksGourp <- merge(allRanksGourp, allRanks[,-c(ncol(allRanks)-1,ncol(allRanks))], by = "row.names")
      rownames(allRanksGourp) <- allRanksGourp[,"Row.names"]
      allRanksGourp <- allRanksGourp[,-1]
  }

  rank.elv <- merge(allRanksGourp, env[,c("Elevation","ForestType")], by = "row.names")
  rank.elv <- rank.elv[order(rank.elv[,"Elevation"]),]
  colnames(rank.elv)[1] <- "Samples"
  
  fname <- paste("max-remained", divs[i], "group-heatmap.pdf", sep = "-")
  plotHeatmapMaxRemDivRank(rank.elv, fname, title=subTitles[i])
  
  rownames(rank.elv) <- rank.elv[,"Samples"]
  rank.elv <- rank.elv[,-c(1,ncol(rank.elv)-1,ncol(rank.elv))]
  printXTable(rank.elv, caption = paste("Ranking sampling plots by removing plots sequentially so as to 
              minimize the loss of overall", divs[i], "diversity among the remaining plots. 
              1 is the most important and removed at the last, 10 is the least important and removed in the beginning."), 
              label = paste0("tab:maxRemained", divs[i], "Group"), file=tableFile)
}



############### in development ##############

q=0
ylab=expression(paste("% of maximum "^"0", D[gamma]))
for (lev in c("gamma", "beta")) {
    
    if (lev == "beta") {
		q=1
		ylab=expression(paste("% of maximum "^"1", D[beta]))
	}
	
	# ranks, last 2 col are mean and sd
	allRanks <- getMaxRemainedDiversityRank(lev, q, taxa.group) 
	# max diversity
	allMaxDiv <- getMaxRemainedDiversity(lev, q, taxa.group)

	if (lev == "gamma") {
	  printXTable(allRanks, caption = "Ranking sampling plots by removing plots sequentially so as to 
			 minimize the loss of overall $\\gamma$ diversity among the remaining plots. 
	     1 is the most important and removed at the last, 10 is the least important and removed in the beginning.", 
	     label = "tab:maxRemainedGamma0", file=tableFile)
	} else {
	  printXTable(allRanks, caption = "Ranking sampling plots by removing plots sequentially so as to 
				maximize the resulting effective $\\beta$ diversity among the remaining plots. 
				1 is the most important and removed at the last, 10 is the least important and removed in the beginning.",
				label = "tab:maxRemainedBeta1", file=tableFile)		
	}
	
	p <- ggplot(gg, aes(x= reorder(variable, value, median, order=TRUE), y=value)) + 
	  geom_boxplot(aes(fill = factor(variable))) + 
	  ggtitle(paste("Cluster", cl)) + xlab("") + ylab("Abundence") + 
	  theme_bw() + guides(fill=FALSE) +
	  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
	        panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) 
	
	######## figure #######
	if(!all(rownames(allRanks) == rownames(allMaxDiv))) 
	  stop(paste("Find different plots ", rownames(allRanks), " != ", rownames(allMaxDiv), sep=""))
	
	pdf(file.path(workingPath, figDir, paste0("max-remained-", lev, q, ".pdf")), width=6, height=5)
	attach(mtcars)
	par(mfrow=c(3,2), oma=c(3,0,0,0))
	
	# 1st col is Rank
	for (expId in 1:(n-1)) {
	  val <- data.frame(x=allRanks[,expId], y=allMaxDiv[,expId]/max(allMaxDiv[,expId]), stringsAsFactors = FALSE)
	  val <- val[order(val$x, decreasing=T),]
	  if (expId == 1) {
	    par(mar=c(4,5,2,1))		
	    plot(val, pch=myshape[expId], col=mypalette[expId], xlim=c(nrow(val),1), ylim=c(0,1), 
	         xlab="number of sites remained", ylab=ylab, main="(a)") 				  		
	  } else {
	    points(val, pch=myshape[expId], col=mypalette[expId]) 
	  }
	  
	  lines(val, lty=2, col=mypalette[expId]) 	
	}  # END for expId
	
	par(usr=c(0,1,0,1),xpd=NA)  
	legend(-1.38, -0.3, ncol=3, legend=matrixNames, pch=as.numeric(myshape), col=mypalette)
	
	invisible(dev.off())  
	

} # END for lev


