# eDNA-pipeline

In development

## Python Scripts
1. run\_PEAR.py is a script to run PEAR paired-end read merger. 

  Zhang et al. 2014, PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics 30:5, 614-620. http://sco.h-its.org/exelixis/web/software/pear/

2. fastq\_primer\_split.py splits fastq files containing multiple amplicons into separate files based on primer sequences. It requires a tab-delimited file with gene names and forward/reverse primer sequences. 

3. relabel\_reads.py adds sample names and other info to sequence ids in fastq files, also to rename fastq files by sample. It requires a tab-delimited file with sequence output IDs and corresponding sample names. 

4. concatenate\_files.py simply concatenates fastq files together.

5. usearch8\_pipeline\_miseq.py runs USEARCH v8 commands for error filtering, OTU clustering, etc. on fastq files.

6. make\_otutable\_from\_uparse.py makes a table of samples by OTU sequence counts based on sequence dereplication data and output of the UPARSE OTU clustering algorithm in USEARCH v8 (used in usearch8\_pipeline\_miseq.py).

7. make\_taxonomy_table_v6.py takes a file of taxonomic paths for each OTU exported from MEGAN metagenome analyser, and outputs a table of taxonomic ranks, down to genus or the otherwise lowest taxonomic rank for each OTU. Taxonomic ranks at order level and above are reorganised according to the unified taxonomy scheme of Ruggiero et al. (2015), resulting in a consistent set of ranks for all taxonomic groups. (This isn't perfect, as a few taxonomic names in NCBI (Cryptomycota, for example) don't occur in Ruggiero et al. (2015), and vice-versa.) Ranks below order level are organised according to the NCBI taxonomy scheme. 

  Ruggiero et al. 2015, A Higher Level Classification of All Living Organisms. http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0119248

## Taxonomy data
1. New\_taxonomy\_from\_PLOSONE_cleaned_up is the unified taxonomy scheme (to order level) that comes from Ruggiero et al. (2015). 

2. NCBI\_taxa\_list.txt is a set of taxonomic names and ranks derived from the NCBI taxonomy scheme. 
