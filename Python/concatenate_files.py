# -*- coding: utf-8 -*-
"""
@author: A Dopheide
"""

import os, shutil, glob, gzip

os.chdir("../Data")
with open("18S_all_sequences.fastq", "wb") as outfile_1:
    for filename in glob.glob('*.fastq.gz'):
        if "18S_" in filename:
            print filename        
            with gzip.open(filename) as readfile:
                shutil.copyfileobj(readfile, outfile_1)
outfile_1.close()
