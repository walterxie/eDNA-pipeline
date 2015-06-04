
* createAllCommunityMatrix.r
  Create community matrix which is the transverse matrix of community matrix in *vegan* package.
  code is tested using output from usearch8.0.1623
  Input: 1) OTUs mapping file "out.up" from USEARCH OTU clustering; 
  Input: 2) OTU representative sequences "otus.fasta" from USEARCH OTU clustering;
  Input: 3) Chimeras sequences chimeras.fasta from Uchime; 
  Input: 4) Duplicate sequences mapping file "derep.uc" from USEARCH de-replication. 
  Input: 5) SraRunTable.txt, a SRA mapping file to map SRA code to subplot name.
  Output: 1) community matrix CSV file, rows are OTUs, columns are samples.  
  Note: USEARCH 8 has a critical error to mix all sample-based data into one pool,
	therefore the sample information of each duplicate read would be lost during 
	de-replication process. In this code, we retrive the sample of each duplicate 
	read from the mapping file derep.uc created by a modified command:
  ```
  $USEARCH -derep_fulllength ./qc/denoised.fasta -fastaout ./qc/derep.fasta -sizeout -uc ./qc/derep.uc
  ```
  to create the community matrix retaining the correct reads' distribution of samples.  

* createAllDiversitiesOTUsTable.r


* createAllRarefactionTable.r


* allDiversitiesOTUs.r


* allElevationAlpha.r


* allElevationDiversitiesByPlots.r


* allGeneCorrolation.r


* allMDSBySubplots.r


* allMaxDiv.r


* allMaxDivHeatmap.r


* allMaxRemainedDiversity.r


* allRarefactions.r


* allRedundancyAnalysis.r


* allSampleCount.r


* allStatistics.r


* allTaxonomyPhylum.r


* allWithinBetweenPlots.r


* avgRanks.r


