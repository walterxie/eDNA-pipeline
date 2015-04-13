#!/usr/bin/env bash
# cd 16S

export USEARCH="/Applications/usearch8"

../scripts/pipelineDerep.sh $1

for i in {90..100}
do
   ../scripts/pipelineOTUsRef.sh $i ../gold.fa
done
