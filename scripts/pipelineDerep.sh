#!/usr/bin/env bash
USEARCH="/Applications/usearch8"

if [ -z "$1" ]; then
    echo "No fastq file name supplied"
    exit 1
fi

if [ ! -d "qc" ]; then
   mkdir qc
fi

# Quality filter, length truncate, covert to FASTA
$USEARCH -fastq_filter $1 -fastq_trunclen 300 -fastq_maxee 3 -fastaout ./qc/reads.fasta

$USEARCH -derep_fulllength ./qc/reads.fasta -fastaout ./qc/derep.fasta -sizeout

# minsize 2 to filter out all size 1  
$USEARCH -sortbysize ./qc/derep.fasta -fastaout ./qc/sorted.fasta -minsize 1
