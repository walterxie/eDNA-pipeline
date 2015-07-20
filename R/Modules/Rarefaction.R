# depend on source("Modules/RarefractionDiversitiesTable.R")

if(!exists("outputRFTable")) stop("output path for rarefraction table is missing !")

raremax <- min(rowSums(communityMatrix))
rareReps <- 10
sampleSizesPoints <- 11 # include 1

sampleSizesSeq <- round(seq(1, raremax, length.out = sampleSizesPoints), digits = 0)

rdiversityTableMean <- NULL
rdiversityTableSe <- NULL
for(sS in sampleSizesSeq){   
   if (is.null(rdiversityTableMean)) {
        rdiversityTableMean <- data.frame(row.names=paste(levels, qs, sep=""), check.names=FALSE)
		rdiversityTableMean <- rdiversityTable(communityMatrix, sS, rareReps)[,1]
		rdiversityTableSe <- data.frame(row.names=paste(levels, qs, sep=""), check.names=FALSE)
		rdiversityTableSe <- rdiversityTable(communityMatrix, sS, rareReps)[,2]
	} else {
		tmpRDT <- rdiversityTable(communityMatrix, sS, rareReps)[,1] 
		rdiversityTableMean <- cbind(rdiversityTableMean, tmpRDT)
		tmpRDT <- rdiversityTable(communityMatrix, sS, rareReps)[,2] 
		rdiversityTableSe <- cbind(rdiversityTableSe, tmpRDT)
	}  
	
	print(paste("Rarefaction take sample size ", sS, ", max = ", sampleSizesSeq[length(sampleSizesSeq)], sep="")) 
}
colnames(rdiversityTableMean) <- paste("size.", sampleSizesSeq, sep="")
colnames(rdiversityTableSe) <- paste("size.", sampleSizesSeq, sep="")

write.csv(rdiversityTableMean, outputRFTable, row.names=FALSE, quote=FALSE)

#LM.rdt <- lm(rdiversityTableMean[1,]~sampleSizesSeq, rdiversityTableMean[1,])

print(paste("write rarefaction table file to ", outputRFTable, sep=""))
