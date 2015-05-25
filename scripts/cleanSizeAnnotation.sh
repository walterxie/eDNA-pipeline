#!/usr/bin/env bash

THRESHOLD=$1

if [ -z "$1" ]; then
    THRESHOLD=97
fi

for dir in "COI-spun" "COI" "ITS" "trnL" "18S" "16S"; do
	cd ./$dir/otus$THRESHOLD
    echo "Clean $(pwd)/$dir.fasta for MEGAN"
#	sed -i.bak 's/;size=\([0-9]*\);/|\1/g' $dir.out
	sed -i.bak 's/;size=\([0-9]*\);/|\1/g' $dir.fasta
	cd ../..
done
