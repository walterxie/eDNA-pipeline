# -*- coding: utf-8 -*-
"""
Python script to split fastq sequence files containing a mixture of amplicons by primers.
Requires a tab-delimited file with gene and forward/reverse primer sequences.
@author: A Dopheide
"""

from Bio import SeqIO, SeqUtils
import os
import gzip # If input files are gzipped

def split_files_by_primer(fastq_file):
    """Splits each fastq file by gene/primer, outputs trimmed sequences to file"""
    with open(fastq_file) as reads_handle:
    #with gzip.open(fastq_file) as reads_handle: # If files are gzip compressed
        for gene, primers in primerlist.items():
            if "_L001_R1_" in str(fastq_file): # For R1 files, only need to remove primer from one end 
                outfilename = fastq_file.split("_L001_R1")[0]
                with open(("%s_%s_R1_trim.fastq") % (outfilename, gene),"w") as output_handle: 
                    found_reads = split_reads(reads_handle, gene, primers)
                    count = SeqIO.write(found_reads, output_handle, "fastq")
                    print("%s: Saved %i %s reads") % (fastq_file, count, gene)
                    logfile.write(("%s \t %s \t %i \n") % (fastq_file, gene, count))
            elif "pear-contigs.assembled" in str(fastq_file): # For merged contig files, need to remove primers from both ends 
                outfilename = fastq_file.split(".assembled")[0]
                with open(("%s_%s_trim.fastq") % (outfilename, gene),"w") as output_handle: 
                    found_reads = split_contigs(reads_handle, gene, primers)
                    count = SeqIO.write(found_reads, output_handle, "fastq")
                    print("%s: Saved %i %s reads") % (fastq_file, count, gene)
                    logfile.write(("%s \t %s \t %i \n") % (fastq_file, gene, count))


def split_reads(reads_handle, gene, primers):
    """Finds primer at one end of sequences only, returns records with primer removed"""
    primer_F = primers[0]
    length_F = len(primer_F) # used to trim primer
    for record in SeqIO.parse(reads_handle, "fastq"):
        primer_search = SeqUtils.nt_search(str(record.seq), primer_F)  # Searches record.seq for DNA sequence (i.e. primer)
        if len(primer_search) > 1 and (primer_search)[1] == 0:  # check if primer found (len > 1) at start of record ([1] == 0)             
            yield record[length_F:]  # Trims primer length from record, returns trimmed record


def split_contigs(reads_handle, gene, primers):
    """Finds primers at both ends of sequences, returns records with both primers removed"""
    primer_F = primers[0]
    length_F = len(primer_F) # used to trim primer
    primer_R = primers[1]
    length_R = len(primer_R) # used to trim primer
    for record in SeqIO.parse(reads_handle, "fastq"):
        primer_search_F = SeqUtils.nt_search(str(record.seq), primer_F)  # Searches record.seq for DNA sequence (i.e. primer)
        #primer_search_R = SeqUtils.nt_search(str(record.seq), primer_R)  # Searches record.seq for DNA sequence (i.e. primer)
        if len(primer_search_F) > 1 and (primer_search_F)[1] == 0:  # Check if primer found (len > 1) at start of record ([1] == 0)
            yield record[length_F:-length_R] # Trims primer lengths from both ends of record, returns trimmed record
#            if gene == "18S":  
#                yield record[primer_length:-21]
#            elif gene == "ShCO1":
#                yield record[primer_length:-26]
#            elif gene == "16S":
#                yield record[primer_length:-16]
#            elif gene == "26S":
#                yield record[primer_length:-19]
#            elif gene == "ITS":
#                yield record[primer_length:-20]


##### Set up filenames etc #####
os.chdir("../Data/")

# Load list of primers (tab-delimited gene names and forward/reverse primer sequences)
primerlist = {}
with open("Barcode_split_primers.txt","r") as primerfile:
    for row in primerfile:
        gene, primer_F, primer_R = row.strip().split("\t")
        primers = (primer_F, primer_R)
        primerlist[gene] = primers
        print ("%s %s") % (gene, primer_F, primer_R)

# Iterate through fastq(gz?) files to be split
listing = os.listdir("../Data/") # Location of files to be split
with open("Primer_split_log.txt","a") as logfile:
    logfile.write("Reads file \t Primer \t Read count \n")
    for file in listing:
        if str(file).endswith(".fastq"):
            print ("Analysing %s" % file)
            split_files_by_primer("%s" % file)
