library(xtable)
library(vegetarian)

diversity_table <- data.frame(row.names=c("gamma", "alpha", "beta"))

diversity_table$'q=0' <- c(
  d(communityMatrix,lev="gamma",q=0),
  d(communityMatrix,lev="alpha",q=0),
  d(communityMatrix,lev="beta",q=0))

diversity_table$'q=1' <- c(
  d(communityMatrix,lev="gamma",q=1),
  d(communityMatrix,lev="alpha",q=1),
  d(communityMatrix,lev="beta",q=1))

diversity_table$'q=2' <- c(
  d(communityMatrix,lev="gamma",q=2),
  d(communityMatrix,lev="alpha",q=2),
  d(communityMatrix,lev="beta", q=2))

colnames(diversity_table) <- c("$q=0$", "$q=1$", "$q=2$")
rownames(diversity_table) <- c("$D_\\gamma(q)$", "$D_\\alpha(q)$", "$D_\\beta(q)$")
#div.xtable <- xtable(diversity_table, caption=paste("Jost Diversities ($q\\in{0,1,2}$) for ",matrixName, ", where species richness q=0, Shannon q=1, Simpson q=2.", sep=""))
#print(div.xtable, sanitize.text.function = function(x){x})


# COMPUTE TURNOVER TABLE
#
#d.turnover <- TurnoverDist(communityMatrix)
#

# COMPUTER HORN-MORISITA OVERLAPS

#library(vegan)
#d.hornMorisita <- vegdist(communityMatrix, method="horn", binary=FALSE)
#d.brayBin <- vegdist(communityMatrix, method="bray", binary=TRUE)


#library(untb)
#cm_counts <- count(colSums(communityMatrix))
#theta <- round(optimal.theta(cm_counts),2)
