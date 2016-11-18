# Hauturu MiSeq pipeline

Two MiSeq runs were generated. The first was a 2 x 250 bp run containing 18S, ITS, and COI-300 amplicons. The second was a 2 x 300 bp run containing 16S, 26S, and COI-650 amplicons. Each sequencing run resulted in fastq files corresponding to R1 and R2 reads for each sample.

1. Preprocessing part A: 
For the 2 x 250 bp run, use PEAR (http://sco.h-its.org/exelixis/web/software/pear/) to overlap the R1 and R2 reads. PEAR doesn't work if any irregularities present in files to be processed (e.g. reads missing due to mismatch with primers), so do this before anything else:
- With all the R1/R2 files in the same folder, use "run\_PEAR.py" to overlap reads with a minimum 50 bp overlap;
- Outputs overlapped sequences as fastq files corresponding to each pair of input R1/R2 files; 

For the 2 x 300 bp run, the reads don't overlap, so just use the (higher-quality) R1 reads. 

2. Preprocessing part B: 
Split the overlapped or R1 reads fastq files by gene and remove primers: 
- use "fastq\_primer\_split\_v3.py";
- Requires a tab-delimited list of gene names and corresponding primer sequences;
- Outputs a fastq file for each sample and gene with sequences trimmed of primers, renamed as Sample-ID\_contigs/reads\_gene\_primertrim.fastq;

Add gene and sample names to sequence ids in each split/trimmed fastq file:
- use "relabel\_reads.py";
- Requires a tab-delimited list of sequencing IDs and corresponding sample IDs;
- Outputs relabelled fastq files, with gene and sample ID added to sequence records (format: |gene\_x|sample\_y); 

Concatenate the trimmed/relabelled fastq files for each gene
- use "concatenate\_files.py" (or just "cat *.fastq > gene\_all\_contigs/reads.fastq");

3. OTU clustering
Filter, trim, dereplicate and cluster OTUs at 97 % identity threshold from the combined fastq file for each gene using Usearch v8.1.1861 for linux64.
- use "usearch8\_pipeline\_miseq\_2015.py" to run usearch commands, summarised below:
1. Trimming and error filtering:
- Use maximum expected error threshold (maxee) of 1; 
- For 18S, ITS and COI-300 overlapped sequences, use minimum lengths of 310, 180 and 300 bp 
respectively:
"usearch -fastq\_filter file -fastq\_minlen 310 -fastq\_maxee 1.0 -fastaout filter\_trim.fasta"
- For 16S, 26S and COI-650 amplicons, truncate at 250 bp, to remove the lowest-quality regions:
"usearch -fastq\_filter file -fastq\_trunclen 250 -fastq\_maxee 1.0 -fastaout filter\_trim.fasta"
2. Dereplicate filtered sequences (use uc output file for mapping reads to OTUs later):
"usearch derep\_fulllength filter\_trim.fasta -fastaout uniques.fasta -sizeout -uc derep.uc" 
3. Sort dereplicated sequences by size for OTU clustering input:
"usearch -sortbysize uniques.fasta -fastaout uniques\_sorted.fasta"
4. Cluster into OTUs at 97 % identity threshold using UPARSE algorithm:
"usearch -cluster\_otus uniques\_sorted.fasta -otus OTUs.fasta -uparseout OTUs\_out.up"
5. Make an OTU table based on uparse OTUs\_out.up and sequence counts from derep.uc files:
- use "make\_otutable\_from\_uparse.py"

4. Get taxonomic identifications for each OTU
BLAST each OTU file against the NCBI nr database. Import the resulting BLAST output files into MEGAN metagenome analyser. From MEGAN, output taxonomic paths for each OTU. Reorganise the taxonomic paths into table of taxonomic ranks based on the all-inclusive taxonomy scheme of Ruggiero et al. 2015, using "make\_taxonomy\_table\_v6.py".   
