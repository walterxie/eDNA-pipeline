# Create community matrix
# code is tested using output from usearch8.0.1623
# Input: 1) OTUs mapping file "out.up" from USEARCH OTU clustering; 
# Input: 2) OTU representative sequences "otus.fasta" from USEARCH OTU clustering;
# Input: 3) Chimeras sequences chimeras.fasta from Uchime; 
# Input: 4) Duplicate sequences mapping file "derep.uc" from USEARCH de-replication. 
# Note: USEARCH 8 has a critical error to mix all sample-based data into one pool,
#       therefore the sample information of each duplicate read would be lost during 
#       de-replication process. In this code, we retrive the sample of each duplicate 
#       read from the mapping file derep.uc created by a modified command:
# $USEARCH -derep_fulllength ./qc/denoised.fasta -fastaout ./qc/derep.fasta -sizeout -uc ./qc/derep.uc
#       to create the community matrix retaining the correct reads' distribution of samples.  
# Input: 5) SraRunTable.txt, a SRA mapping file to map SRA code to subplot name.
# Output: 1) community matrix CSV file, rows are OTUs, columns are samples.  
# Note: communityMatrix here is the transverse matrix of community matrix in vegan package.
# 
# Author: Walter Xie
# Accessed on 22 May 2015

#source("http://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
#browseVignettes("ShortRead")

library(ShortRead)

# change config below
sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
setwd(sourcePath)

workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
experiments <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name   
matrixNames <-  experiments # for better naming

n <- length(experiments) 

# hard code for SRA data format
# get column (sample) index given a label and column names of CM
# label has no size annotation
whichSampleIn <- function(samplesCM, label) {
    sampleN <- gsub("(SRR[0-9]*).[0-9]*", "\\1", label)
    sampleIdOfCM <- which(samplesCM %in% sampleN) # which(colnames(communityMatrix) %in% "SRR1706107")
    return(sampleIdOfCM)
}

writeCM <- function(communityMatrix, experiment, otuThre) {
	if(sum(rowSums(communityMatrix)==0) | sum(colSums(communityMatrix)==0)) # use which(rowSums(communityMatrix)==0) to check
		print( paste("Warning : community matrix contains zero-read column/row !") )

	print(paste("Create community matrix for ", matrixNames[expId], " (", otuThre, "%) : rows = ", nrow(communityMatrix), 
		", columns = ", ncol(communityMatrix), ", total reads = ", sum(communityMatrix), 
		", singleton = ", sum(rowSums(communityMatrix)==1), ", coupleton = ", sum(rowSums(communityMatrix)==2),
		"; min sample per plot = ", min(colSums(communityMatrix)), "; max = ", max(colSums(communityMatrix)), sep=""))
	
	outputCM <- paste(workingPath, experiment, "/otus", otuThre, "/", experiment, "-", otuThre, ".csv", sep="")
	write.csv(communityMatrix, outputCM, quote=FALSE)
}

# main
sra_table <- read.table(paste(workingPath, "data/SraRunTable.txt", sep=""), sep="\t", header=T, row.names=1, check.names=FALSE)
sra_table[,1] <- gsub("CM30C30", "Plot9", sra_table[,1], ignore.case = T)
sra_table[,1] <- gsub("LB1", "Plot10", sra_table[,1], ignore.case = T)
sra_table <- sra_table[order(rownames(sra_table)),]

for (expId in 1:n) {	
	# UC format http://drive5.com/usearch/manual/opt_uc.html
	derep.uc <- read.table(paste(workingPath, experiments[expId], "/qc/derep.uc", sep=""), sep="\t")
	derep.uc<-derep.uc[derep.uc[,1]=="H",] # only keep duplicates id
	
	# recommend to replace loop by otuThre=97
	for (otuThre in 90:100) {
		###### final OTUs after removing chimeras from Uchime ######
		otus.fasta <- readFasta(paste(workingPath, experiments[expId], "/otus", otuThre, sep=""), "^otus.fasta$")
		chimeras.fasta <- readFasta(paste(workingPath, experiments[expId], "/otus", otuThre, sep=""), "^chimeras.fasta$")

		if (length(chimeras.fasta) <= 0) 
			stop(paste("Chimeras file does not load properly ! ", length(chimeras.fasta), sep=""))

		final.fasta <- otus.fasta[!(id(otus.fasta) %in% id(chimeras.fasta)),]
		
		if (length(otus.fasta) <= length(final.fasta)) 
			stop(paste("Chimeras are not removed properly ! ", length(otus.fasta), " vs. ", length(final.fasta), sep=""))

		#TODO how to rename the label
		#final.fasta <- gsub(";size=([0-9]*);", "|\\1", final.fasta)		

		writeFasta(final.fasta, paste(workingPath, experiments[expId], "/otus", otuThre, "/", experiments[expId], ".fasta", sep=""))
		        
		###### retrive rows (OTUs) & columns (samples) ######
        # all OTUs 
		otus <- id(final.fasta)
		# remove size annotation
		otus <- gsub(";size=([0-9]*);", "", otus)
		# extract samples, e.g.
		samples <- sort(unique(gsub("(SRR[0-9]*).[0-9]*", "\\1", otus)))
		
		print(paste(matrixNames[expId], " (", otuThre, "%) : OTUs = ", length(otus), ", samples = ", length(samples)))
		
        ###### count total chimeras ######   
        # UP format http://drive5.com/usearch/manual/opt_uparseout.html
		# due to strange error for up file, use fill = T, header = F	
		out.up <- read.table(paste(workingPath, experiments[expId], "/otus", otuThre, "/out.up", sep=""), sep="\t", fill = T, header = F)
		
		chimeras.up.len <- length(out.up[out.up[,2]=="chimera", 1])
		otus.up.len <- length(out.up[out.up[,2]=="otu", 1])
		members.up.len <- length(out.up[out.up[,2]=="match", 1])
		
		print( paste("Total unique reads =", (members.up.len+otus.up.len) , ", during OTU clustering.") )
		print( paste("Total chimeras =", (chimeras.up.len+length(chimeras.fasta)) , ", where", chimeras.up.len, "from OTU clustering and", 
					length(chimeras.fasta), "from Uchime.") )
        
        if (length(otus) != (otus.up.len - length(chimeras.fasta))) 
			print( paste("Warning : up file has incorrect number of OTUs after removing chimeras :", (otus.up.len - length(chimeras.fasta))) )
		
		###### create CM ######        		
		communityMatrix <- matrix(0,nrow=length(otus),ncol=length(samples))
		rownames(communityMatrix) <- otus
		colnames(communityMatrix) <- samples
		
		###### fill in CM ######
		pb <- txtProgressBar(min=1, max=nrow(communityMatrix), style = 3)
		for (i in 1:nrow(communityMatrix)) {
			# txtProgressBar
			setTxtProgressBar(pb, i)

			# OTU representative sequence label
			otuN <- rownames(communityMatrix)[i]		
			otuSampleIdOfCM <- whichSampleIn(colnames(communityMatrix), otuN)		
			communityMatrix[i,otuSampleIdOfCM] <- communityMatrix[i,otuSampleIdOfCM] + 1
			# duplicate reads for OTU representative sequence
			duplicatesN <- derep.uc[derep.uc[,10]==otuN,9] 
			for (dupN in duplicatesN) {
				duplicateSampleIdOfCM <- whichSampleIn(colnames(communityMatrix), dupN)	
				communityMatrix[i,duplicateSampleIdOfCM] <- communityMatrix[i,duplicateSampleIdOfCM] + 1
			} # end for dupN

			# unique reads belonging to this OTU
		    membersN <- out.up[(out.up[,2]=="match" & out.up[,5]==otuN), 1]
		    # remove size annotation
			membersN <- gsub(";size=([0-9]*);", "", membersN)
		    for (memN in membersN) {
				memberSampleIdOfCM <- whichSampleIn(colnames(communityMatrix), memN)	
				communityMatrix[i,memberSampleIdOfCM] <- communityMatrix[i,memberSampleIdOfCM] + 1
				
				# duplicate reads for OTU member sequence
				duplicatesN <- derep.uc[derep.uc[,10]==memN,9] 
				for (dupN in duplicatesN) {
					duplicateSampleIdOfCM <- whichSampleIn(colnames(communityMatrix), dupN)	
					communityMatrix[i,duplicateSampleIdOfCM] <- communityMatrix[i,duplicateSampleIdOfCM] + 1
				} # end for dupN
			} # end for memN
		} # end for i
		close(pb)
		
		communityMatrix <- communityMatrix[,order(colnames(communityMatrix))]
		selectedSamples <- sra_table[which(rownames(sra_table) %in% colnames(communityMatrix)),]
		print("Selected samples from SRA mapping file : ")
		print(selectedSamples)
			
		gene <- unique(selectedSamples[,2])
		if (gene != experiments[expId]) 
			stop(paste("Sample name is not matched to a correct data set ! ", gene))
		
		if (nrow(selectedSamples) != ncol(communityMatrix)) 
			stop(paste("Sample name in community matrix not exist in SRA mapping table ! "))
		
		print( paste("Replace sample name in community matrix [", paste(colnames(communityMatrix),collapse=", "), 
			"] by [", paste(selectedSamples[,1],collapse=", "), "].") )
		
		colnames(communityMatrix) <- selectedSamples[,1]
		
		writeCM(communityMatrix, experiments[expId], otuThre)		
	} # end for otuThre

} # end for expId
