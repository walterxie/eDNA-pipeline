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
  
  merge.df <- merge(phylo.rare.df, env[,c("Plot","Elevation","ForestType")], by = "row.names", sort = FALSE)
  # rename Row.names
  colnames(merge.df)[1] <- "Samples"
  
  melt.df <- melt(merge.df, id=c("Samples","Plot","Elevation","ForestType"))
  melt.df <- melt.df[complete.cases(melt.df),]
  # remove prefix size.
  melt.df$variable  <- gsub("^.*?\\.","",melt.df$variable) 
  
  melt.df$variable <- as.numeric(melt.df$variable)
  melt.df$Samples = factor(melt.df$Samples, unique(melt.df$Samples))
  melt.df$ForestType = factor(melt.df$ForestType, sort(unique(melt.df$ForestType)))
  
  p <- ggplot(data=melt.df, aes(x=variable, y=value, group=Samples, shape=ForestType, colour=Elevation)) + 
    scale_shape_manual(values=0:nlevels(melt.df$ForestType)) +
    xlab("Rarefaction Size") + ylab("Mean of Rooted PD") + 
    ggtitle(paste(matrixNames[expId], taxa.group)) +     
    geom_line(size=.5) +     
    geom_point(size=3, fill="white") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank()) 
  
  fname <- paste("phylorare", matrixNames[expId], postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-")
  pdf(file.path(workingPath, figDir, paste(fname, "pdf", sep = ".")), width=9, height=6)
  print(p)
  invisible(dev.off())   
}
