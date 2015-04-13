#!/usr/bin/env bash

for dir in "COI-spun" "COI" "ITS" "trnL" "16S" "18S"; do
	if [ ! -d "$dir" ]; then
	    echo "Cannot find folder $dir"
        exit 1
	fi
	
	if [ ! -f CMCreator0.1.jar ]; then
		echo "Jar file CMCreator not found!"
		exit 1
	fi

	if [ ! -f SraRunTable.txt ]; then
		echo "Mapping file SraRunTable.txt not found!"
		exit 1
	fi

	for i in {90..100}; do
	   	if [ ! -d "$dir/otus$i" ]; then
			echo "Cannot find folder $dir/otus$i"
			exit 1
		fi
        
        if [ ! -f $dir/otus$i/out.up ]; then
			echo "UP file out.up not found!"
			exit 1
		fi

		if [ ! -f $dir/otus$i/chimeras.fasta ]; then
			echo "Chimeras file chimeras.fasta not found!"
			exit 1
		fi

        java -jar CMCreator0.1.jar -working $dir/otus$i -out $dir-$i.csv -chi chimeras.fasta -sra data/SraRunTable.txt -overwrite out.up
	done
	
done
