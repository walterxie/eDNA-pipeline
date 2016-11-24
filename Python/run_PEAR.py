# -*- coding: utf-8 -*-
"""

Python script to run PEAR 
http://sco.h-its.org/exelixis/web/software/pear/
Zhang et al. 2014. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics 30:5, 614-20.  

@author: A Dopheide
"""

import os, subprocess

listing = os.listdir("../Data/") # Location of R1/R2 data files

def find_matching_file (R1_file):
    R1_file_start = R1_file[:9]
    print R1_file
    for fastq_file in listing:
        if fastq_file.startswith(R1_file_start) and "R2_001" in fastq_file:
            return fastq_file

def run_PEAR(R1_file, R2_file, outfile):
    start_PEAR = "pear -f %s -r %s -v 50 -o %s" % (R1_file, R2_file, outfile) # Run PEAR, require 50 bp overlap
    print starting_PEAR
    subprocess.call(start_PEAR.split(), stdout = log, stderr = subprocess.STDOUT)

with(open("PEAR_log.txt", "a") as log:
    for f in listing:
        if str(f).endswith("R1_001.fastq"):
            print "R1 file: %s" % str(f)
            R1_file = f
            R2_file = find_matching_file(R1_file)
            print "R2 file: %s" % R2_file
            sample_id = infile[:9]
            outfile = "%s_pear-contigs_" % sample_id
            print "Outfile %s" % outfile
            run_PEAR(R1_file, R2_file, outfile)