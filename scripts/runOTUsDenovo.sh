#!/usr/bin/env bash
# cd 18S

export USEARCH="/Applications/usearch8"

../scripts/pipelineDerep.sh $1

for i in {90..100}
do
   ../scripts/pipelineOTUsDenovo.sh $i
done
