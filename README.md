# eDNA-pipeline

Evaluating a multigene environmental DNA approach for comprehensive biodiversity assessment

## Dataset

1. BioProject in NCBI:
  http://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA267737

2. 20 BioSamples (including lat-long, elevation, temperature):
  http://www.ncbi.nlm.nih.gov/biosample?Db=biosample&DbFrom=bioproject&Cmd=Link&LinkName=bioproject_biosample&LinkReadableName=BioSample&ordinalpos=1&IdsFromResult=267737

3. Soil environmental sequences (454) in SRA:
  http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=SRP050103

4. Elevations:
  plot_elevations.txt

5. Invertebrate (leaf litter + pitfall traps):
  CO1_Invertebrate_Pitfall_1526_OTUs.csv, CO1_Leaf_Litter_1526_OTUs.csv

6. Vegetation survey data:
  Hauturu (Little Barrier island) in Jan 2011 (http://nvs.landcareresearch.co.nz).
  diameters_pilot_species_1_2_0.csv, saplings_pilot_species_1_2_0.csv, seedlings_pilot_species_1_2_0.csv

7. Bird counts:
  birds_pilot_species_1_2_0.csv

8. SraRunTable.txt : 
  a SRA mapping file to map SRA code to subplot name.
  
9. Soil chemistry:
  LJ12027.txt

## Folder structure in working path 

1. A working folder, such as ./pipeline

2. Folders for each data set (genes), such as ./pipeline/16S

3. Folders for deconvolution, such as ./pipeline/16S/deconvoluted

4. Folders for quality control, such as ./pipeline/16S/qc

5. Folders for each OTU threshold, such as ./pipeline/16S/otus97

For example,
```
pipeline
 |___ 16S
 |     |___ deconvoluted
 |     |___ qc
 |     |___ otu?? (e.g. 97)
 |___ 18S
 ...  
```


## Download 454 sequences 

1. Download soil environmental sequences (454) from SRA:
  http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=SRP050103

2. Convert sra format into fastq using SRA Toolkit, such as:
  ```
  /sratoolkit.2.4.3-mac64/bin/fastq-dump SRR1720812.sra
  ```
  SRR1720812.fastq will be ready in the same folder.

3. Combine all fastq files from the same marker into one file:
  ```
  cat *.fastq > 16S.fastq
  ```
  Alternatively, use scripts *downloadData.sh* and *prepareData.sh* as discribed below.


## UPARSE 

1. Download USEARCH (http://www.drive5.com/usearch/download.html)

2. Setup USEARCH by copying it to /Applications, and create a link:
  ```
  ln -s usearch8.0.??? usearch8
  ```

3. Create a working folder, such as ./pipeline

4. Copy all scripts in to working folder, such as pipeline/scripts

5. Download data
  ```
  scripts/downloadData.sh 
  scripts/prepareData.sh
  ```

6. We strongly recommend to use error corrector for 454 data. In this pipeline, 
we choose Acacia (http://www.nature.com/nmeth/journal/v9/n5/abs/nmeth.1990.html). 
Please install Acacia and change its path in the script *pipelineDerep.sh* before go to step 7.
Note that you need to give enough memory according to the size of the largest cluster from your data. 
We recommend to assign at least 50 GB memory to Acacia (-Xmx50g) for our 16S and 10 GB for the rest of datasets. 

7. Go into each gene folder to run script for either denovo chimera filtering or reference filtering (16S only), e.g.
  ```
  cd 16S
  ../scripts/runOTUsRef.sh ./deconvoluted/16S.fastq 
  ```
  or 
  ```
  cd COI
  ../scripts/runOTUsDenovo.sh ./deconvoluted/COI.fastq 
  ```
  Note: the reference dataset gold.fa (http://www.drive5.com/usearch/manual/cmd_uchime_ref.html) is required by *runOTUsRef.sh*.


## Generate Community Matrix 

Note: USEARCH 8 mixes all sample-based data into one pool,
therefore the sample information of each duplicate read would be lost during 
de-replication process. In *createAllCommunityMatrix.r*, we retrive the sample of each duplicate 
read from the mapping file derep.uc created by a modified command:
```
$USEARCH -derep_fulllength ./qc/denoised.fasta -fastaout ./qc/derep.fasta -sizeout -uc ./qc/derep.uc
```
to create the community matrix retaining the correct reads' distribution of samples.   


##  BLAST 

Use *cleanSizeAnnotation.sh* to clean up the size annotation in the sequence label created by USEARCH, 
which will cause MEGAN input error. 


## Community Matrix Analysis 

1. Change source path for pipeline and working path for data into your local path.

2. Run *createAllCommunityMatrix.r* first to create community matrices.

3. Run *createAllDiversitiesOTUsTable.r* second to get the rarefaction table and OTU threshold table. This is time-consuming.
If you do not need clustering through different thresholds, you could use faster script *createAllRarefactionTable.r* to generate 
the rarefaction table at 97% threshold only. 

4. Run *all???.r* one by one to get all figures and tables.
