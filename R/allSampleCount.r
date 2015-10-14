library(ggplot2)

# change config below
#figDir <- "figures"
#sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
#setwd(sourcePath)
#workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
#matrixNames <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name   

if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("isPlot")) isPlot = FALSE # by subplot
if(!exists("otuThr")) otuThr = 97
if(!exists("levels")) levels = rep(c("gamma","alpha","beta"),3)
if(!exists("qs")) qs = rep(0:2,each=3)

n <- length(matrixNames) 
mypalette <- rainbow(n)
myshape <- seq(0, (0 + n-1))

individualName <- "read"
sampleName <- "site"
speciesName <- "OTU"
speciesNamePlural <- "OTUs"

#source("Modules/init.R", local=TRUE)

######## Reads OTUs #######
labls = c("reads rest", "150 most abundant OTUs", "OTUs 1 read", "OTUs 2 reads", "OTUs >= 3")
cat=c(rep("reads",2), rep("OTUs",3))
myPalette <- c("#fdbb84", "#e34a33", "#e0f3db", "#a8ddb5", "#43a2ca")

cat("\nGraph: reads counts percentage bar chart including singletons: isPlot =", isPlot, ", otuThr =", otuThr, "\n") 
# always plot singleton
for (expId in 1:n) {	
	communityMatrix <- getCommunityMatrixT(expId, isPlot, FALSE)

	source("Modules/SampleCounts.R", local=TRUE)

	if (expId == 1) {
		readsOTUs <- data.frame(cat=cat, labls=labls,
			value=c((individual_count-threshold_individual_count) / individual_count * 100, threshold_individual_count / individual_count * 100, 
			length(cs[cs == 1]) / length(cs) * 100, length(cs[cs == 2]) / length(cs) * 100, length(cs[cs != 1 & cs != 2]) / length(cs) * 100))
		readsOTUs$region <- matrixNames[expId]
	
	} else {
		readsOTUs_tmp <- data.frame(cat=cat, labls=labls,
			value=c((individual_count-threshold_individual_count) / individual_count * 100, threshold_individual_count / individual_count * 100, 
			length(cs[cs == 1]) / length(cs) * 100, length(cs[cs == 2]) / length(cs) * 100, length(cs[cs != 1 & cs != 2]) / length(cs) * 100))
		readsOTUs_tmp$region <- matrixNames[expId]
		readsOTUs <- rbind(readsOTUs, readsOTUs_tmp) 
	}

}

readsOTUs$region = factor(readsOTUs$region,matrixNames)
readsOTUs$labls = factor(readsOTUs$labls,labls)

pdf(paste(workingPath, figDir, "/", postfix("reads-counts", isPlot, FALSE, sep="-"), ".pdf", sep = ""), width=8, height=5)	

print( ggplot(readsOTUs, aes(x = region, y = value, fill = labls)) + geom_bar(stat="identity", position='stack') +
	 theme_bw() + facet_grid( ~ cat) + ylab("Percentage") + scale_fill_manual(breaks=rev(labls), values= myPalette) + 
	 theme(legend.position="top", legend.title=element_blank(), axis.title.x=element_blank()) )

invisible(dev.off()) 


######## Sample counts #######

subTitles <- c("(a)","(b)","(c)","(d)","(e)","(f)")

cat("\nGraph: sample counts bar chart: rmSingleton =", rmSingleton, ", isPlot =", isPlot, ", otuThr =", otuThr, "\n") 

pdf(paste(workingPath, figDir, "/", postfix("sample-counts", isPlot, rmSingleton, sep="-"), ".pdf", sep = ""), width=6, height=9)	
attach(mtcars)
par(mfrow=c(3,2))	

for (expId in 1:n) {	
  # isPlot determines to use which matrix file in init
  communityMatrix <- getCommunityMatrixT(expId, isPlot, rmSingleton)
    
  source("Modules/SampleCounts.R", local=TRUE)
    
	if (expId > 4) {		
		xlab=paste("Number of ",sampleName,"s crossed", sep="")
		par(mar=c(4,4,3,2)) 
	} else {
		xlab=""
		par(mar=c(4,4,3,2))
	}
	
	if (expId %% 2 != 0) {		
		ylab=paste("Number of ",speciesNamePlural, "", sep="")		
	} else {
		ylab=""
	}
			
  barplot(log10(table(sampleCounts)), xlab=xlab, ylab=ylab, main=paste(subTitles[expId], matrixNames[expId], sep = " "), yaxt="n")   
    
  aty <- axTicks(2)
	ylabels <- sapply(aty,function(i) 10^i)
	#labels <- sapply(aty,function(i) as.expression(bquote(2^ .(i)))  )
	axis(2,at=aty,labels=ylabels)    
}
invisible(dev.off()) 


