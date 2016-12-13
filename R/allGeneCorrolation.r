# Figure 5: Mantel test (lower triangle) and Procrustes test (upper triangle) correlations for comparisons

# corrs <- getMantelAndProcrustes(input.names)
getMantelAndProcrustes <- function(input.names, metric="jaccard",
                    genes.taxa=list(list("16S","all"),list("18S","all"),list("26S","all"),
                                    list("ITS","all"),list("ShCO1","all"),list("FolCO1","all")),
                    order.by=c("16S all", "18S all", "26S all", "ITS all", "COI-300 all", "COI-650 all") ) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)
  
  cm.by.subplot.list <- getCommunityList(genes=input.names, genes.taxa=genes.taxa, by.plot=F, drop.taxa=TRUE )
  cat("\n")
  
  dissim <- ComMA::getDissimilarityList(cm.by.subplot.list, metric=metric)
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

# corrs$procrustes
plotProcrustes.allOTUs <- function(procrustes, env.subplot) {
  cat("Max elevation is ", max(env.subplot[,"Elevation"]), "\n")
  theme_set(theme_bw(base_size=8))
  # Figure S8
  pS8.list <- ComMA::plotProcrustes(procrustes$proc.list, env.subplot, proc.list.pairs=procrustes$pairs, 
                                    colour.id="Elevation", limits=c(0,650))
  pS8 <- ComMA::grid_arrange_shared_legend(pS8.list[[1]], input.list=T, ncol=3, nrow=5, 
                                           legend.position="right", widths=c(0.8, 0.15, 0))
}

# corrs2$procrustes
plotAllProcrustes <- function(procrustes, env.subplot) {
  theme_set(theme_bw(base_size=8))
  # Figure 7
  p7.list <- ComMA::sublistByPairs(procrustes$proc.list, procrustes$pairs, 
                                   subset.pairs=list(list("16S bacteria","18S protists"),list("16S bacteria","18S fungi"),list("16S bacteria","18S animals"),
                                                     list("16S bacteria","26S fungi"),list("16S bacteria","COI-300 animals"),list("18S protists","18S fungi"),
                                                     list("18S protists","18S animals"),list("18S protists","26S fungi"),list("18S protists","COI-300 animals"),
                                                     list("18S fungi","18S animals"),list("18S fungi","26S fungi"),list("18S fungi","COI-300 animals"),
                                                     list("18S animals","26S fungi"),list("18S animals","COI-300 animals"),list("26S fungi","COI-300 animals")) )
  p7 <- ComMA::plotProcrustes(p7.list$sub.list, env.subplot, proc.list.pairs=p7.list$sub.pairs, 
                              colour.id="Elevation", limits=c(0,650))
  # Figure S9
  pS9.list <- ComMA::sublistByPairs(procrustes$proc.list, procrustes$pairs, 
                                    subset.pairs=list(list("18S fungi","26S fungi"),list("18S fungi","ITS fungi"),list("18S fungi","COI-300 fungi"),
                                                      list("26S fungi","ITS fungi"),list("26S fungi","COI-300 fungi"),list("ITS fungi","COI-300 fungi")) )
  pS9 <- ComMA::plotProcrustes(pS9.list$sub.list, env.subplot, proc.list.pairs=pS9.list$sub.pairs, 
                               colour.id="Elevation", limits=c(0,650))
  # Figure S10
  S10.list <- ComMA::sublistByPairs(procrustes$proc.list, procrustes$pairs, 
                                    subset.pairs=list(list("18S protists","26S protists"),list("18S protists","COI-300 protists"),list("18S protists","COI-650 protists"),
                                                      list("26S protists","COI-300 protists"),list("26S protists","COI-650 protists"),list("COI-300 protists","COI-650 protists")) )
  S10 <- ComMA::plotProcrustes(S10.list$sub.list, env.subplot, proc.list.pairs=S10.list$sub.pairs, 
                               colour.id="Elevation", limits=c(0,650))
  # Figure S11
  pS11.list <- ComMA::sublistByPairs(procrustes$proc.list, procrustes$pairs, 
                                    subset.pairs=list(list("18S animals","26S animals"),list("18S animals","COI-300 animals"),list("18S animals","COI-650 animals"),
                                                      list("26S animals","COI-300 animals"),list("26S animals","COI-650 animals"),list("COI-300 animals","COI-650 animals")) )
  pS11 <- ComMA::plotProcrustes(pS11.list$sub.list, env.subplot, proc.list.pairs=pS11.list$sub.pairs, 
                                colour.id="Elevation", limits=c(0,650))
  # Figure S12
  pS12.list <- ComMA::sublistByPairs(procrustes$proc.list, procrustes$pairs, 
                                    subset.pairs=list(list("18S fungi","18S protists"),list("18S fungi","18S animals"),list("18S protists","18S animals")) )
  pS12 <- ComMA::plotProcrustes(pS12.list$sub.list, env.subplot, proc.list.pairs=pS12.list$sub.pairs, 
                                colour.id="Elevation", limits=c(0,650))
  # Figure S13
  pS13.list <- ComMA::sublistByPairs(procrustes$proc.list, procrustes$pairs, 
                                    subset.pairs=list(list("26S fungi","26S protists"),list("26S fungi","26S animals"),list("26S protists","26S animals")) )
  pS13 <- ComMA::plotProcrustes(pS13.list$sub.list, env.subplot, proc.list.pairs=pS13.list$sub.pairs, 
                                colour.id="Elevation", limits=c(0,650))
  # Figure S14
  pS14.list <- ComMA::sublistByPairs(procrustes$proc.list, procrustes$pairs, 
                                    subset.pairs=list(list("COI-300 fungi","COI-300 protists"),list("COI-300 fungi","COI-300 animals"),list("COI-300 protists","COI-300 animals")) )
  pS14 <- ComMA::plotProcrustes(pS14.list$sub.list, env.subplot, proc.list.pairs=pS14.list$sub.pairs, 
                                colour.id="Elevation", limits=c(0,650))
  # Figure S15
  pS15.list <- ComMA::sublistByPairs(procrustes$proc.list, procrustes$pairs, 
                                    subset.pairs=list(list("COI-650 protists","COI-650 animals")) )
  pS15 <- ComMA::plotProcrustes(pS15.list$sub.list, env.subplot, proc.list.pairs=pS15.list$sub.pairs, 
                                colour.id="Elevation", limits=c(0,650))
  
  # return a list of gtable
  list(gt7=ComMA::grid_arrange_shared_legend(p7[[1]], input.list=T, ncol=3, nrow=5, legend.position="right", widths=c(0.8, 0.15, 0)),
       gtS9=ComMA::grid_arrange_shared_legend(pS9[[1]], input.list=T, ncol=3, nrow=2, legend.position="right", widths=c(0.8, 0.15, 0)),
       gtS10=ComMA::grid_arrange_shared_legend(S10[[1]], input.list=T, ncol=3, nrow=2, legend.position="right", widths=c(0.8, 0.15, 0)),
       gtS11=ComMA::grid_arrange_shared_legend(pS11[[1]], input.list=T, ncol=3, nrow=2, legend.position="right", widths=c(0.8, 0.15, 0)),
       gtS12=ComMA::grid_arrange_shared_legend(pS12[[1]], input.list=T, ncol=3, nrow=1, legend.position="right", widths=c(0.8, 0.15, 0)),
       gtS13=ComMA::grid_arrange_shared_legend(pS13[[1]], input.list=T, ncol=3, nrow=1, legend.position="right", widths=c(0.8, 0.15, 0)),
       gtS14=ComMA::grid_arrange_shared_legend(pS14[[1]], input.list=T, ncol=3, nrow=1, legend.position="right", widths=c(0.8, 0.15, 0)),
       gtS15=ComMA::grid_arrange_shared_legend(pS15[[1]][[1]], ncol=2, nrow=1, legend.position="right", widths=c(1,0.3))
  )
}
