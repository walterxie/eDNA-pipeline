# -*- coding: utf-8 -*-
"""
Python script to relabel fastq files and sequences by samples.
Requires a tab-delimited table of well ids and corresponding sample names.
@author: A Dopheide
"""

from Bio import SeqIO
import os

os.chdir("../Data")

# Sample ID file is a tab-delimited table of NZGL wells/IDs and corresponding sample names
sample_ids = {}
with open("NZGL01401_run2_18S_ITS_sample_ID.txt", "r") as id_file: 
    for row in id_file:
        sample, name = row.strip().split("\t")
        sample_ids[sample] = name
        print ("%s %s") % (sample, name)

# ShCO1 done separately to 18S/ITS due to differing sample configuration (three different wells)
sample_ids_ShCO1 = {}
with open("NZGL01401_run2_ShCO1_sample_ID.txt", "r") as id_file_ShCO1:
    for row in id_file_ShCO1:
        sample, name = row.strip().split("\t")
        sample_ids_ShCO1[sample] = name
        print ("%s %s") % (sample, name)

# Load a (trimmed) fastq file, get sample name, relabel reads with sample name, 
# output relabelled reads to new file:
listing = os.listdir("../Data/")
for infile in listing:
    if str(infile).endswith("trim.fastq")
        print ("Relabelling %s" % infile)
        filestart = infile[:6] # get NZGL well number
        if ("18S") in str(infile) or ("ITS") in str(infile):
            sample_name = sample_ids.get(filestart)  # Get matching sample name
        elif("ShCO1") in str(infile):
            sample_name = sample_ids_ShCO1.get(filestart)  # Get (different) matching sample name for ShCO1
        print sample_name
        name_parts = str.split(infile, "_")  # Get parts for relabelling
        with open(infile) as reads_handle:
            if "R1" in str(infile):
                output_handle = open(("%s_%s_%s_%s.fastq") % (name_parts[1], sample_name, name_parts[2], "R1-reads"), "w")  # for read files (not contigs)
            elif "contigs" in str(infile):
                output_handle = open(("%s_%s_%s_%s.fastq") % (name_parts[1], sample_name, name_parts[2], "pear-contigs"), "w")  # for contig files (not reads)
            for record in SeqIO.parse(reads_handle, "fastq"):
                record.id = (("%s|gene_%s|sample_%s") % (record.id, name_parts[2], sample_name))  # add locus and sample name to record label
                SeqIO.write(record, output_handle, "fastq")
        output_handle.close()
        #reads_handle.close()
