#!/usr/bin/env bash

dataset=( "COI-spun" "16S" "18S" "COI" "ITS" "trnL" )
#SRR1706131
accession=( "SRR1720793" "SRR1706032" )
runs=20

# COI-spun SRR index increases 1 by 1
i=0
if [ ! -d "${dataset[$i]}" ]; then
   mkdir ${dataset[$i]}
fi

cd ./${dataset[$i]}
if [ ! -d "data" ]; then
   mkdir ./data
fi

cd ./data
echo "${dataset[$i]}"
for (( a = 0 ; a < $runs ; a++ )) do
	acc=$((${accession[0]:3:10}+$a))
	newAcc="SRR$acc"
	echo "download $newAcc"
	ftp ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${newAcc:0:6}/$newAcc/$newAcc.sra
done
cd ..
cd ..

# every 4 SRR indexes belong to a same mark
for (( i = 1 ; i < ${#dataset[@]} ; i++ )) do
	if [ ! -d "${dataset[$i]}" ]; then
	   mkdir ${dataset[$i]}
	fi

	cd ./${dataset[$i]}
	if [ ! -d "data" ]; then
	   mkdir ./data
	fi
	
	cd ./data
	echo "${dataset[$i]}"
	for (( a = 0 ; a < $runs ; a++ )) do
		acc=$((${accession[1]:3:10}+$a*5+i-1))
		newAcc="SRR$acc"
		echo "download $newAcc"
		ftp ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${newAcc:0:6}/${newAcc}/${newAcc}.sra
	done
	cd ..
	cd ..
done