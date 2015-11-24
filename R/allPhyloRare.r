# phylorare: http://davidnipperess.blogspot.co.nz/search/label/Software


if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("tableFile")) stop("table file is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("taxa.group")) taxa.group="assigned"
if(!exists("isPlot")) isPlot = FALSE # by subplot

n <- length(matrixNames) 

cat("\nAnalysis: plot the rarefied Phylogenetic Diversity of eDNA community. \n")

env <- getSampleMetaData(isPlot)
env[,"ForestType"] <- gsub(":.*", "", env[,"ForestType"], ignore.case = T)
env[,"ForestType"] <- gsub("x", "unknown", env[,"ForestType"], ignore.case = T)

for (expId in 1:(n-1)) {
  ### eDNA ###
  phylo.rare.df <- getPhyloRareTable(expId, isPlot, rmSingleton, taxa.group)
  
  if (is.null(phylo.rare.df)) {
    cat("\nSkip", taxa.group, "subset from", matrixNames[expId], "because of no tree file.\n")
    next
  }
  
  merge.df <- merge(phylo.rare.df, env[,c("Elevation","ForestType")], by = "row.names", sort = FALSE)
  # rename Row.names
  colnames(merge.df)[1] <- "Samples"
  
  melt.df <- melt(merge.df, id=c("Samples","Elevation","ForestType"))
  melt.df <- melt.df[complete.cases(melt.df),]
  # remove prefix size.
  melt.df$variable  <- gsub("^.*?\\.","",melt.df$variable) 
  
  melt.df$variable <- as.numeric(melt.df$variable)
  melt.df$Samples = factor(melt.df$Samples, unique(melt.df$Samples))
  #sort(unique(melt.df$ForestType))
  melt.df$ForestType = factor(melt.df$ForestType, c("VS2","VS3","VS5","WF7","WF9","WF11","WF12","WF13","MF20","unknown"))
  
  end.lines <- aggregate(cbind(variable, value) ~ Samples, melt.df, max)
  end.lines <- merge(end.lines, unique(melt.df[,c("Samples","Elevation","ForestType")]), by = "Samples", sort = FALSE)
  
  p <- ggplot(data=melt.df, aes(x=variable, y=value, group=Samples, colour=Elevation)) + 
    scale_shape_manual(values=c(15,16,17,0,1,2,5,6,3,4)) +
    scale_colour_gradientn(colours = c("blue", "orange")) +
    scale_x_sqrt(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) + 
    #scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), label=scientific_10) +
    xlab("Rarefaction Size") + ylab("Mean of Rooted PD") + 
    ggtitle(paste(matrixNames[expId], taxa.group)) +     
    geom_line(size=.5) +     
    geom_point(data=end.lines, aes(shape=ForestType, colour=Elevation), size=3) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank()) 
  
  fname <- paste("phylorare", matrixNames[expId], postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-")
  pdf(file.path(workingPath, figDir, paste(fname, "pdf", sep = ".")), width=9, height=6)
  print(p)
  invisible(dev.off())   
}
