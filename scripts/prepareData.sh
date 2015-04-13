#!/usr/bin/env bash

dataset=( "COI-spun" "16S" "18S" "COI" "ITS" "trnL" )
fastqdump="/Applications/sratoolkit.2.4.5-mac64/bin/fastq-dump"

for ds in "${dataset[@]}"; do
	if [ ! -d "$ds" ]; then
		echo "cannot find dataset $ds !"
		exit 0
	fi
	cd ./$ds
	
	if [ ! -d "data" ]; then
		echo "cannot find data folder in $ds !"
		exit 0
	fi
	if [ ! -d "deconvoluted" ]; then
		mkdir ./deconvoluted
	fi
		
	cd ./data
	$fastqdump *.sra
	cd ..
	
	echo "creating ${ds}.fastq ..."
	cat ./data/*.fastq > ./deconvoluted/${ds}.fastq
	
	cd ..
done