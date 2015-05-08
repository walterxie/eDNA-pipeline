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
  trees_pilot_species_1_2_0.csv, seedlings_pilot_species_1_2_0.csv

7. Bird counts:
  birds_pilot_species_1_2_0.csv


## Folder structure in working path 

1. A working folder, such as *./pipeline*

2. Folders for each data set (genes), such as *./pipeline/16S*

3. Folders for deconvolution, such as *./pipeline/16S/deconvoluted*

4. Folders for quality control, such as *./pipeline/16S/qc*

5. Folders for each OTU threshold, such as *./pipeline/16S/otus97*


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

1. Download USEARCH http://www.drive5.com/usearch/download.html

2. Setup USEARCH by copying it to /Applications, and create a link:
  ```
  ln -s usearch8.0.??? usearch8
  ```

3. Create a working folder, such as *./pipeline*

4. Copy all scripts in to working folder, such as *pipeline/scripts*

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
  Note: the reference dataset gold.fa (http://www.drive5.com/usearch/manual/cmd_uchime_ref.html) is required.


## Generate Community Matrix 

1. Create a data folder under the working folder, such as *./pipeline/data*. And download modified SraRunTable.txt from 

2. Stay in the working folder, such as ./pipeline
  ```
  java -jar CMCreator0.1.jar -working COI/otus97 -out COI-97.csv -chi chimeras.fasta -sra data/SraRunTable.txt -overwrite out.up
  ```


## Community Matrix Analysis 

1. Change source path for pipeline and working path for data into your local path.

2. Run createAllDiversitiesOTUsTable.r first to get the rarefraction table and OTU threshold table. This is time-consuming.

3. Run all*.r one by one to get all figures and tables.
