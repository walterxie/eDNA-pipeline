#!/usr/bin/env bash

THRESHOLD=$1

if [ -z "$1" ]; then
    THRESHOLD=97
fi

for dir in "COI-spun" "COI" "ITS" "trnL" "16S" "18S"; do
	cd ./$dir/otus$THRESHOLD
    echo "Clean $(pwd)/otus.fasta for MEGAN"
#	sed -i.bak 's/;size=\([0-9]*\);/|\1/g' otus.out
	sed -i.bak 's/;size=\([0-9]*\);/|\1/g' otus.fasta
	cd ../..
done
