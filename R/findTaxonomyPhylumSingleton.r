# change config below
sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/"
setwd(sourcePath)
workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
matrixNames <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") 

if(!exists("matrixNames")) stop("matrix names are missing !")

rmSingleton <- FALSE
n <- length(matrixNames) 
otuThr = 97

source("Modules/init.R", local=TRUE)

allTaxaSingleton <- NULL
for (expId in 1:n) {
	communityMatrix <- init(expId, otuThr, "-by-subplot")
	singletons <- which(colSums(communityMatrix)==1)

	taxaToReadsName <- read.table(paste(workingPath, "data/", matrixNames[expId], "-ex.txt", sep=""), sep="\t", stringsAsFactors = FALSE)
	colnames(taxaToReadsName) <- c("OTU","taxa")
	# SRR1706072.75769|1	Actinobacteria
	tipLabels <- strsplit(taxaToReadsName[,1], "\\|")
	taxaToReadsName[,1]<-unlist(tipLabels)[2*(1:length(tipLabels))-1]

	matched <- match(colnames(communityMatrix)[singletons], taxaToReadsName[,1])
	matched <- matched[!is.na(matched)]

	taxaAssgSingleton <- aggregate(.~taxa, data=taxaToReadsName[matched,], FUN=length)
	
	if (expId == 1) {
		allTaxaSingleton <- taxaAssgSingleton
	} else {
		allTaxaSingleton <- merge(allTaxaSingleton, taxaAssgSingleton, by = "taxa", all = TRUE)
	}
	colnames(allTaxaSingleton)[expId+1] <- matrixNames[expId]	
}

allTaxaSingleton[is.na(allTaxaSingleton)] <- 0

allTaxaSingleton$taxa <- gsub(" <phylum>", "", allTaxaSingleton$taxa, ignore.case = T)
allTaxaSingleton$taxa <- gsub(" \\(.*\\)", "", allTaxaSingleton$taxa, ignore.case = T)
allTaxaSingleton$taxa <- gsub("Nucleariidae and Fonticula group", "Nucleariidae & Fonticula", allTaxaSingleton$taxa, ignore.case = T)
allTaxaSingleton$taxa <- gsub("environmental samples <Rhizaria>", "unclassified Rhizaria", allTaxaSingleton$taxa, ignore.case = T)

uclBa1 <- grep("candidate division", allTaxaSingleton$taxa, ignore.case = FALSE)
uclBa2 <- grep("Bacteria", allTaxaSingleton$taxa, ignore.case = FALSE)

uclBa <- c("unclassified Bacteria", colSums(allTaxaSingleton[c(uclBa1,uclBa2),2:7]))

allTaxaSingleton <- allTaxaSingleton[-c(uclBa1,uclBa2),]
allTaxaSingleton <- rbind(allTaxaSingleton, uclBa)
allTaxaSingleton <- allTaxaSingleton[order(as.numeric(rownames(allTaxaSingleton))),]

write.table(allTaxaSingleton, paste(workingPath, "data/taxonomy97singleton.txt", sep=""), sep="\t", row.names = FALSE)
