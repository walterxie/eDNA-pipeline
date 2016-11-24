# -*- coding: utf-8 -*-
""""
Pipeline for processing Illumina Miseq output fastq files using USEARCH V8.  

@author: Andrew Dopheide

"""

import os
import subprocess

def filter_derep_reads(input_file, label, trunclen, maxee, log):
    """Filter and dereplicate reads"""
    # 1. Trim and quality filter reads
    print ("Filtering %s\nlabel: %s") % (input_file, label)
    filter_trim = "%s -fastq_filter %s -fastq_trunclen %i -fastq_maxee %.1f \
        -fastaout %s_filter_trim.fasta" % (usearchpath, input_file, trunclen, maxee, label)
    #log.write('\nFilter and trim reads:\n' + str(filter_trim) + '\n')
    subprocess.call(filter_trim.split(), stdout=log, stderr=subprocess.STDOUT)
    # 2. Dereplicate reads
    derep = "%s -derep_fulllength %s_filter_trim.fasta \
        -fastaout %s_uniques.fasta -sizeout -uc %s_derep.uc" % (usearchpath, label, label, label)
    #log.write('\nDereplicate filtered/trimmed reads:\n' + str(derep) + '\n')
    subprocess.call(derep.split(), stdout=log, stderr=subprocess.STDOUT)


def filter_derep_contigs(input_file, label, minlen, maxee, log):    
    """Filter and dereplicate contigs"""
    # 1. Trim and quality filter contigs
    print ("Filtering %s\nlabel: %s") % (input_file, label)
    filter_trim = "%s -fastq_filter %s -fastq_minlen %i -fastq_maxee %.1f \
        -fastaout %s_filter_trim.fasta" % (usearchpath, input_file, minlen, maxee, label)
    #log.write('\nFilter and trim reads:\n' + str(filter_trim) + '\n')
    subprocess.call(filter_trim.split(), stdout=log, stderr=subprocess.STDOUT)
    # 2. Dereplicate contigs
    derep = "%s -derep_fulllength %s_filter_trim.fasta \
        -fastaout %s_uniques.fasta -sizeout -uc %s_derep.uc" % (usearchpath, label, label, label)
    #log.write('\nDereplicate filtered/trimmed contigs:\n' + str(derep) + '\n')
    subprocess.call(derep.split(), stdout=log, stderr=subprocess.STDOUT)


def get_OTUs(label, log):
    """Sort by size and cluster OTUs, make OTU table from Usearch8 derep.uc and uparse outfiles"""
    # 3: Sort by size
    print ("sorting by size")
    sortbysize = "%s -sortbysize %s_uniques.fasta \
        -fastaout %s_uniques_sorted.fasta" % (usearchpath, label, label)
    #log.write('\nSort reads by size, minsize 1:\n' + str(sortbysize) + '\n')
    subprocess.call(sortbysize.split(), stdout=log, stderr=subprocess.STDOUT)
    # 4: Cluster OTUs
    print ("clustering otus")
    clusterOTUs = "%s -cluster_otus %s_uniques_sorted.fasta -otus %s_OTUs.fasta \
        -uparseout %s_OTUs_out.up -sizein -sizeout" % (usearchpath, label, label, label)
    #log.write('\nCluster reads, including singletons, into OTUs at 97%:\n' + str(clusterOTUs) + '\n')
    subprocess.call(clusterOTUs.split(), stdout=log, stderr=subprocess.STDOUT)
    # 5: Make OTU table
    #makeOTUtable = "python make_otutable_from_uparse.py %s" % (label)
    #subprocess.call(makeOTUtable.split(), stdout=log, stderr=subprocess.STDOUT)


##### Main code starts here #####
usearchpath = "../usearch8.0.1623_i86linux64"
sys.path.append('../python/')
from make_otutable_from_uparse import make_otutable

maxee = 1  # Error trimming threshold
trunclen = 250  # Length trimming cutoff for reads

listing = os.listdir("../seqdata/")
for input_file in listing:
    labelpart = str.split(input_file, ".fastq")[0]
    label = ("%s_ME%.1f" % (labelpart, maxee))
    log = open(("%s_usearch_log.txt" % label), "a")
    if str(input_file).endswith("contigs.fastq"):
        if "18S" in str(input_file):
            minlen = 310
        elif "ITS" in str(input_file):
            minlen = 180
        elif "ShCO1" in str(input_file):
            minlen = 300
        filter_derep_contigs(input_file, label, minlen, maxee, log)
    elif str(input_file).endswith("reads.fastq"):
        filter_derep_reads(input_file, label, trunclen, maxee, log)
    get_OTUs(label, log)
    up_file = '%s_OTUs_out.up' % label
    derep_file = '%s_derep.uc' % label
    make_otutable(label, up_file, derep_file)