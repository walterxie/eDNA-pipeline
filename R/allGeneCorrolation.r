


# mantel <- getMantel()
getMantelAndProcrustes <- function(input.names, metric="jaccard",
                    genes.taxa=list(list("16S","all"),list("18S","all"),list("26S","all"),
                                    list("ITS","all"),list("ShCO1","all"),list("FolCO1","all")) ) {
  if (missing(input.names)) 
    source("R/init.R", local=TRUE)
  
  cm.list <- getCommunityList(genes=input.names, genes.taxa=genes.taxa, by.plot=F, 
                              col.ranks=c("superkingdom", "kingdom"), drop.taxa=TRUE )
  cat("\n")
  
  dissim <- ComMA::getDissimilarityList(cm.list, metric=metric)
  mantel <- ComMA::mantelComparison(dissim$dist.list)
  procrustes <- ComMA::procrustesComparison(dissim$dist.list)
  
  mantel.tri <- ComMA::getTriMatrix(mantel$m.df) # Mantel stats
  prot.tri <- ComMA::getTriMatrix(procrustes$prot.df) # Procrustes stats
  corrs <- ComMA::combineTriMatrix(mantel.tri, prot.tri)
  
  list( corrs=corrs, mantel=mantel, procrustes=procrustes, dist.list=dissim$dist.list, 
        metric=metric, msg="Mantel in upper triangle, Procrustes in lower"  )
}


# m.mantel <- matrix(0,nrow=(n+m),ncol=(n+m))
# colnames(m.mantel) <- c(matrixNames, matrixNamesNo454)
# rownames(m.mantel) <- c(matrixNames, matrixNamesNo454)
# m.signif <- matrix(0,nrow=(n+m),ncol=(n+m))
# colnames(m.signif) <- c(matrixNames, matrixNamesNo454)
# rownames(m.signif) <- c(matrixNames, matrixNamesNo454)
# 
# for (expId1 in 1:(n+m-1)) {
# 	for (expId2 in (expId1+1):(n+m)) {	     
# 		if (expId1 <= n) {
# 			d.comm1 <- beta1DisList[[ expId1 ]]
# 			exp.comm1 <- matrixNames[expId1]
# 		} else {
# 			d.comm1 <- beta1No454DisList[[ expId1-n ]]
# 			exp.comm1 <- matrixNamesNo454[expId1-n]
# 		}
# 
# 		if (expId2 <= n) {
# 			d.comm2 <- beta1DisList[[ expId2 ]]
# 			exp.comm2 <- matrixNames[expId2]
# 		} else {
# 			d.comm2 <- beta1No454DisList[[ expId2-n ]]
# 			exp.comm2 <- matrixNamesNo454[expId2-n]
# 		}
# 
# 		# hard code for missing data in invertebrates
# 		if (length(d.comm1) > length(d.comm2)) {
# 			m.comm <- as.matrix(d.comm1) 
# 			d.comm1 <- as.dist(m.comm[c(-7,-8),c(-7,-8)]) 
# 		} else if (length(d.comm1) < length(d.comm2)) {
# 			m.comm <- as.matrix(d.comm2) 
# 			d.comm2 <- as.dist(m.comm[c(-7,-8),c(-7,-8)]) 	     
# 		}
# 
# 		####### mantel test #######
# 		mantel.comm <- mantel(d.comm1, d.comm2, permutations=4999) 
# 
# 		m.mantel[expId2, expId1] <- mantel.comm$statistic
# 		m.mantel[expId1, expId2] <- mantel.comm$statistic
# 		m.signif[expId2, expId1] <- mantel.comm$signif
# 		m.signif[expId1, expId2] <- mantel.comm$signif
# 
# 		print(paste(exp.comm1, " vs. ", exp.comm2, ", mantel statistic r = ", round(mantel.comm$statistic, 3), "; significance = ", mantel.comm$signif, "; permutations = ", mantel.comm$permutations, sep=""))
# 	}
# }
# 
# ####### MDS of correlation #####
# # Classical MDS
# for (expId1 in 1:(n+m)) {
#     m.mantel[expId1, expId1] <- 1
# }
# d.mantel <- 1-m.mantel # euclidean distances between the rows
# 
# fit <- cmdscale(d.mantel, eig=TRUE, k=2) # k is the number of dim
# 
# # plot solution 
# x <- fit$points[,1]
# y <- fit$points[,2]
# 
# labels = gsub("invertebrates", "inverts", row.names(fit$points), ignore.case = T)
# 
# pdf(paste(workingPath, figDir, "/mds-pairewise-cm-", otuThr, ".pdf", sep = ""), width=5, height=5)
# plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS (1-Correlation)",	type="n")
# text(x, y, labels = labels, cex=.7, xpd=TRUE)
# invisible(dev.off()) 
# 
# output <- paste(workingPath, "mds-pairewise-mantel-", otuThr, ".csv", sep = "")
# # tsv cannot display the 1st cell of columns
# write.csv(d.mantel, output, quote=FALSE)
# 
# # no birds
# m.mantel1 <- m.mantel[-(n+m),-(n+m)] 
# d.mantel <- 1-m.mantel1 # euclidean distances between the rows
# fit <- cmdscale(d.mantel,eig=TRUE, k=2) # k is the number of dim
# 
# # plot solution 
# x <- fit$points[,1]
# y <- fit$points[,2]
# 
# labels = gsub("invertebrates", "inverts", row.names(fit$points), ignore.case = T)
# 
# pdf(paste(workingPath, figDir, "/mds-no-birds-pairewise-cm-", otuThr, ".pdf", sep = ""), width=5, height=5)
# plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS (1-Correlation)", type="n")
# text(x, y, labels = labels, cex=.7, xpd=TRUE)
# invisible(dev.off()) 
# 
# ####### mantel test #######
# df.mantel <- data.frame(m.mantel, check.names=FALSE)
# for (expId1 in 1:(n+m-1)) {
# 	for (expId2 in (expId1+1):(n+m)) {	     
#         df.mantel[expId1, expId2] <- c("")
#         tmp <- paste(round(m.mantel[expId2, expId1], 3), " (", m.signif[expId2, expId1], ")", sep="")
#         df.mantel[expId2, expId1] <- c(tmp)
# 	}
# }
# for (expId1 in 1:(n+m)) {
#     df.mantel[expId1, expId1] <- c("")
# }
# df.mantel <- df.mantel[-1,-ncol(df.mantel)]
#  
# print(xtable(df.mantel, caption = "Pairwise community matrix correlations of effective $\\beta$ diversity within and between 
# 	the eDNA  data sets and traditional data sets, Mantel statistic $r$, and their significance in parentheses 
# 	using Mantel's test based on 4999 permutations.", label = "tab:geneCorr", caption.placement = "top"), 
# 	sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)

