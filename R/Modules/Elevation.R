#### by elevation ####
# TODO: cannot figure out how to lookup value easily in R
# hard code sample_count == 8 to deal with CO1 Invertebrate Pitfall CM lacking of plot 7 & 8
# elevPlot <- elevPlot[-7] twice

# validate input
if(!exists("inputElevation")) stop("inputElevation is missing !")

elev <- read.table(inputElevation, header=TRUE)

print(paste("sample_count = ", sample_count, sep=""))

if (sample_count == 10) {
  elevPlot <- elev$Elevation.metres[seq(1,20, by=2)]
} else if (sample_count == 8) {
  elevPlot <- elev$Elevation.metres[seq(1,20, by=2)]
  elevPlot <- elevPlot[-7]
  elevPlot <- elevPlot[-7]
} else {
  elevPlot <- elev$Elevation.metres[1:20]
}

elevPlotDist <- dist(elevPlot)

d.brayBin <- vegdist(communityMatrix, method="bray", binary=TRUE)
mantel.brayBin <- mantel(elevPlotDist, d.brayBin, permutations=4999)

pdf(paste("figures/elevation-brayBin-", matrixName, ".pdf", sep = ""))
plot(as.vector(elevPlotDist), as.vector(d.brayBin), ylim=c(0,1), pch=16, col="red", 
  xlab="elevation difference (metres)", ylab="Whittaker's beta")
lm.d.brayBin <- lm(as.vector(d.brayBin)~as.vector(elevPlotDist))
abline(lm.d.brayBin)
invisible(dev.off()) 

d.hornMorisita <- vegdist(communityMatrix, method="horn", binary=FALSE)
mantel.hornMorisita <- mantel(elevPlotDist, d.hornMorisita, permutations=4999)

pdf(paste("figures/elevation-hornMorisita-", matrixName, ".pdf", sep = ""))
plot(as.vector(elevPlotDist), as.vector(d.hornMorisita), ylim=c(0,1), pch=16, 
  col="red", xlab="elevation difference (metres)", ylab="Horn-Morisita overlap")
lm.d.hornMorisita <- lm(as.vector(d.hornMorisita)~as.vector(elevPlotDist))
abline(lm.d.hornMorisita)
invisible(dev.off()) 

#turnover_table <- TurnoverDist(communityMatrix)
#mantel(elevPlotDist, turnover_table, permutations=4999)
#
#plot(as.vector(elevPlotDist), as.vector(turnover_table), ylim=c(0,1))
#lm.tunover_table <- lm(as.vector(turnover_table)~as.vector(elevPlotDist))
#abline(lm.tunover_table)