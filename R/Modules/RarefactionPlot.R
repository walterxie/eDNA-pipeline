
# pathFileStem = filePath + fileStem
plotRarefactionTable <- function(rfT.df, env, pathFileStem, ylab="Mean of Rooted PD", verbose=T) {
  merge.df <- merge(rfT.df, env[,c("Elevation","ForestType")], by = "row.names", sort = FALSE)
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
  end.lines$ForestType = factor(end.lines$ForestType, c("VS2","VS3","VS5","WF7","WF9","WF11","WF12","WF13","MF20","unknown"))
  
  p <- ggplot(data=melt.df, aes(x=variable, y=value, group=Samples, colour=Elevation)) + 
    scale_shape_manual(values=c(15,16,17,0,1,2,5,6,3,4)) +
    scale_colour_gradientn(colours = c("blue", "orange")) +
    scale_x_sqrt(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) + 
    #scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), label=scientific_10) +
    xlab("Rarefaction Size") + ylab(ylab) + 
    ggtitle(paste(matrixNames[expId], taxa.group)) +     
    geom_line(size=.5, alpha=0.75) +     
    #geom_point(data=end.lines, aes(shape=ForestType, colour=Elevation), size=3) +
    geom_text(data=end.lines, aes(label=Samples), size=3, alpha=0.5, hjust=-0.1, vjust=-0.6) + #, position=position_jitter(width=-1, height=3)) + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank()) 
  
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  
  pdf(paste(pathFileStem, "pdf", sep = "."), width=8, height=6)
  #print(p)
  print(grid.draw(gt))
  invisible(dev.off())   
}
