
* downloadData.sh		
Download soil environmental sequences (454) from SRA

* prepareData.sh
Combine to one fastq file

* runOTUsDenovo.sh
UPARSE pipeline using denovo chimeras filtering  

* runOTUsRef.sh
UPARSE pipeline using reference chimeras filtering, only applicable to 16S.  

* pipelineDerep.sh	
UPARSE pipeline dereplication, sorting part

* pipelineOTUsDenovo.sh	
UPARSE pipeline OTU clustering and denovo chimeras filtering part

* pipelineOTUsRef.sh
UPARSE pipeline OTU clustering and reference chimeras filtering part, only applicable to 16S.

* createCM.sh		
Create community matrix from UP file.

* cleanSizeAnnotation.sh	
Clean ";size" annotation from label because it will cause a problem for MEGAN.
