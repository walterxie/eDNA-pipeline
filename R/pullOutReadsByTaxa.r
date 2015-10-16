# Create a figure: the number of reads assigned to a pre-defined rank level and also high-level taxa, and coloured by another given rank level.
# Input: 1) taxonomy table from MEGAN, 1st column is OTUs, the last column is taxonomy, the middles are samples (plots/subplots);
# Input: 2) community matrix, but it is the transverse matrix of community matrix in vegan package;
# Input: 3) rankLevel, the number of reads assigned to;
# Input: 4) groupLevel, for colouring, and must higher than rankLevel;
# Output: 1) bar chart.
# 
# Author: Walter Xie
# Accessed on 16 June 2015

#source("http://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
#browseVignettes("ShortRead")

library(ShortRead)

# change config below
workingPath <- "~/svn/compevol/research/NZGenomicObservatory/MiSeq/DOCt1/"
matrixNames <-  c("16S", "18S", "26S", "ITS", "FolCO1", "ShCO1") # only for cm file name and folder name   
taxaFiles <- c("DOCx9_16S_R1_reads_OTUs_num_nt_taxatable.txt", "DOCx9_18S_PEARcontigs_ME1.0_OTUs_num_nt_taxatable.txt", 
				"DOCx9_26S_R1_reads_OTUs_num_nt_taxatable.txt", "DOCx9_ITS_PEARcontigs_ME1.0_OTUs_num_nt_taxatable.txt", 
				"DOCx9_FolCO1_R1_reads_OTUs_num_nt_taxatable.txt", "DOCx9_ShCO1_PEARcontigs_ME1.0_OTUs_num_nt_taxatable.txt")

otuFiles <- c("16S-R1_DOC_reads_OTUs_num.fasta", "DOCx9_18S_PEARcontigs_ME1.0_OTUs_num.fasta", 
				"26S-R1_DOC_reads_OTUs_num.fasta", "DOCx9_ITS_PEARcontigs_ME1.0_OTUs_num.fasta", 
				"FolCO1-R1_DOC_reads_OTUs_num.fasta", "DOCx9_ShCO1_PEARcontigs_ME1.0_OTUs_num.fasta")

#belongTo = "Insecta" 
belongToList = c("Viridiplantae", "Insecta") 

minOTUs = 100

for (belongId in 1:length(belongToList)) {
	belongTo = belongToList[belongId]
	cat(paste("\nChoose", belongTo, ", set minOTUs =", minOTUs, " :\n"))

	for (expId in 1:length(matrixNames)) {
		cat(paste("\nLoad data for", matrixNames[expId], " :\n"))
	
		##### load taxa path #####
		inputTaxa <- paste(workingPath, "Taxonomy_tables/", taxaFiles[expId], sep="")
		taxaPaths <- read.table(inputTaxa, header=TRUE, row.names=1, sep = "\t", stringsAsFactors=FALSE)  

		print( paste("Load", nrow(taxaPaths) , "OTUs in total from taxa path file.") )

		##### filter out rows not belong to given taxa belongTo ##### 
		taxaPaths <- taxaPaths[which(grepl(belongTo, taxaPaths[,1])),] # taxaPaths[,1] is taxa path separated by ;
	
#		if (belongTo == "Metazoa") { # non Arthropoda
#			taxaPaths <- taxaPaths[-which(grepl("Arthropoda", taxaPaths[,1])),]
#		}

		print( paste("Select", nrow(taxaPaths) , "OTUs from taxa path file.") )

		if (nrow(taxaPaths) > minOTUs) {
			##### load otus #####
			otus.fasta <- readFasta(paste(workingPath, "OTUs/", otuFiles[expId], sep="") )

			print( paste("Load", length(otus.fasta) , "OTUs in total from OTUs file.") )

			# cannot rename sequence label
			rownames(taxaPaths) <- gsub(";size_", ";size=", rownames(taxaPaths))
			final.fasta <- otus.fasta[match( rownames(taxaPaths), as.character(id(otus.fasta)) ), ]

			print( paste("Pull out", length(final.fasta) , "OTU representative reads matching selected OTUs.") )

			##### write selected sequences #####
			final.fasta.file <- paste(workingPath, "data/", belongTo, "-", matrixNames[expId], "-otus.fasta", sep="")
			writeFasta(final.fasta, final.fasta.file)
			
			#### replace all ";size=" to ";size_" for MEGAN problem #####
			final <- readLines(final.fasta.file)
			repl <- gsub( ";size=", ";size_", final )
			cat(repl, file=final.fasta.file, sep="\n")
  
		} else {
			print( "Not enough OTUs, jump to the next of the loop" )
		}

	} #END for expId
} #END for belongId
