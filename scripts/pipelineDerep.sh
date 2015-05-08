#!/usr/bin/env bash
USEARCH="/Applications/usearch8"
ACACIA="java -Xmx10g -jar ../../Acacia-1.53/acacia-01.53.b0.jar"

if [ -z "$1" ]; then
    echo "No fastq file name supplied"
    exit 1
fi

if [ ! -d "qc" ]; then
   mkdir qc
fi

# Quality filter, length truncate, covert to FASTA
$USEARCH -fastq_filter $1 -fastq_trunclen 300 -fastq_maxee 3 -fastqout ./qc/postqc.fastq

# 454 Error corrector, trim length to 300 here
$ACACIA -c denoise.config
# rename
cp ./qc/denoise_all_tags.seqOut ./qc/denoised.fasta

$USEARCH -derep_fulllength ./qc/denoised.fasta -fastaout ./qc/derep.fasta -sizeout

# minsize 2 to filter out all size 1  
$USEARCH -sortbysize ./qc/derep.fasta -fastaout ./qc/sorted.fasta -minsize 1
