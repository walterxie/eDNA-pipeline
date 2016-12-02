# Figure 5: Mantel test (lower triangle) and Procrustes test (upper triangle) correlations for comparisons

# corrs <- getMantelAndProcrustes(input.names)
getMantelAndProcrustes <- function(input.names, metric="jaccard",
                    genes.taxa=list(list("16S","all"),list("18S","all"),list("26S","all"),
                                    list("ITS","all"),list("ShCO1","all"),list("FolCO1","all")),
                    order.by=c("16S all", "18S all", "26S all", "ITS all", "COI-300 all", "COI-650 all") ) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)
  
  cm.list <- getCommunityList(genes=input.names, genes.taxa=genes.taxa, by.plot=F, 
                              col.ranks=c("superkingdom", "kingdom"), drop.taxa=TRUE )
  cat("\n")
  
  dissim <- ComMA::getDissimilarityList(cm.list, metric=metric)
  mantel <- ComMA::mantelComparison(dissim$dist.list)
  procrustes <- ComMA::procrustesComparison(dissim$dist.list)
  
  mantel.tri <- ComMA::getTriMatrix(mantel$m.df, order.by=order.by) # Mantel stats
  prot.tri <- ComMA::getTriMatrix(procrustes$prot.df, order.by=order.by) # Procrustes stats
  corrs <- ComMA::combineTriMatrix(mantel.tri, prot.tri)

  rownames(corrs) <- gsub(" all", "", rownames(corrs))
  colnames(corrs) <- gsub(" all", "", colnames(corrs))
  rownames(mantel.tri) <- gsub(" all", "", rownames(mantel.tri))
  colnames(mantel.tri) <- gsub(" all", "", colnames(mantel.tri))
  rownames(prot.tri) <- gsub(" all", "", rownames(prot.tri))
  colnames(prot.tri) <- gsub(" all", "", colnames(prot.tri))
  
  mantel.s.tri <- ComMA::getTriMatrix(mantel$m.df, value="sign", order.by=order.by) # Mantel stats
  prot.s.tri <- ComMA::getTriMatrix(procrustes$prot.df, value="sign", order.by=order.by) # Procrustes stats
  signs <- ComMA::combineTriMatrix(mantel.s.tri, prot.s.tri)
  
  list( corrs=corrs, signs=signs, mantel.tri=mantel.tri, prot.tri=prot.tri,
        mantel=mantel, procrustes=procrustes, dist.list=dissim$dist.list, 
        metric=metric, msg="Mantel test (lower triangle) and Procrustes test (upper triangle) correlations"  )
}

# plotMantelAndProcrustes(corrs)
plotMantelAndProcrustes <- function(corrs, text.repel=T, gene.levels=c("16S", "18S", "26S", "ITS", "COI-300", "COI-650")) {
  min.cor <- min(corrs$corrs[corrs$corrs>0])
  cat("min correlation is ", min.cor, "max is ", max(corrs$corrs), "\n.")

  m.p.tri <- corrs$corrs
  m.p.tri$gene <- rownames(m.p.tri)
  m.p.tri[m.p.tri==0] <- NA
  
  p.hm <- ComMA::ggHeatmap(m.p.tri, melt.id="gene", title="", legend.title="Correlations", 
                        x.levels=gene.levels, y.levels=rev(gene.levels), label.digits=3,
                        high = "#f46d43", mid = "#ffffbf", low = "#3288bd", 
                        midpoint = mean(c(0.4, 1)), limit = c(0.4, 1), 
                        breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) )
  
  p.m.mds <- ComMA::ggNMDSPlot(dist(corrs$mantel.tri), text.or.point=1, text.size=5, text.repel=text.repel,
                               title="Mantel test", title.add.stress=F, title.hjust=0.5) + 
    expand_limits(x = c(-1, 1), y=c(-1, 1))
  
  p.p.mds <- ComMA::ggNMDSPlot(dist(corrs$prot.tri), text.or.point=1, text.size=5, text.repel=text.repel,
                               title="Procrustes test", title.add.stress=F, title.hjust=0.5) + 
    expand_limits(x = c(-1, 1), y=c(-1, 1))
  
  list( heatmap=p.hm, mantel.mds=p.m.mds, prot.mds=p.p.mds, 
        msg="Mantel test (lower triangle) and Procrustes test (upper triangle) correlations"  )
}

printMantelAndProcrustes <- function(corrs, label = "tab:gene:comp", file.xtable=NULL, invalid.char=FALSE,
                                     caption = paste("Mantel and Procrustes test statistics for comparisons", 
                                                     "between the overall OTU community structure detected in", 
                                                     "56 soil samples according to meta-barcoding analysis of different amplicons.")) {
  df <- ComMA::prettyNumbers(corrs$corrs, digits = 3)
  df[df==0] <- ""
  align.v <- rep("r", ncol(df) + 1)
  ComMA::printXTable(df, align = align.v, label = label, file = file.xtable, 
                     invalid.char=invalid.char, caption = caption)
}

