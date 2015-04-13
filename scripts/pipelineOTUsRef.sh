#!/usr/bin/env bash
#USEARCH="/Applications/usearch8"

if [ -z "$1" ]; then
    echo "No OTU threshold supplied"
    exit 1
fi

if [ -z "$2" ]; then
    echo "No chimeras database supplied"
    exit 1
fi

THRESHOLD=$1
PCT=$((100-$THRESHOLD))

REF_DB=$2

if [ ! -d "otus$THRESHOLD" ]; then
   mkdir otus$THRESHOLD
fi

$USEARCH -cluster_otus ./qc/sorted.fasta -otu_radius_pct $PCT -otus ./otus$THRESHOLD/otus.fasta -sizein -sizeout -uparseout ./otus$THRESHOLD/out.up
# -log ./otus$THRESHOLD/OTU$THRESHOLD.log

# only for 16S
$USEARCH -uchime_ref ./otus$THRESHOLD/otus.fasta -db $REF_DB -strand plus -uchimeout ./otus$THRESHOLD/results.uchime -chimeras ./otus$THRESHOLD/chimeras.fasta
#-nonchimeras ./otus$THRESHOLD/otusNoChi.fasta 
