# -*- coding: utf-8 -*-
"""
Aceesed on 19/11/2015

@author: Andrew D

Modified by Walter Xie
"""
import os
#import re
from os.path import expanduser

#filepath = ('H:/My Documents/PhD Research PFR folder/')
filepath = ('%s/Projects/FishGutMicrobiomes/data/' % expanduser('~'))
#filepath = ('C:/Documents and Settings/Andrew/Desktop/')

#os.chdir("%sWaitakere_combined_stuff" % filepath)
os.chdir('%s16s-blast/' % filepath)
#os.chdir('%sLBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/Taxa_tables/' % filepath)

# Organise taxonomic reference data
tax_ref = {}
l = []
tax_ref_file = open('%s/WorkSpace/eDNA-pipeline/data/New_taxonomy_from_PLOSONE_2015_fixed.txt' % expanduser('~'), 'r')
n = 1
for row in tax_ref_file:
    print n
    if n is 1:
        n += 1
    else:
        fields = row.strip('\n').split('\t')
        last = fields[13] # Making dict excludes duplicates
        print last
        tax_ref[last] = row
        l.append(last)
        n += 1

listing = os.listdir('%s16s-blast/' % filepath)
#listing = os.listdir('%sLBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/Taxa_tables/' % filepath)

def get_match(t, tax_ref):
    m = [value for key, value in tax_ref.items() if t.lower() == key.lower()]
    #print len(m)
    print "m: "
    print m
    if len(m) is 1:  # i.e. exact match
        match = m[0]
    elif len(m) > 1:
        match = m[0]
    elif len(m) is 0:
        m2 = [value for key, value in tax_ref.items() if t.lower() in key.lower()]            
        if len(m2) is 1:
            match = m2[0]
        elif len(m2) > 1:
            match = m2[0]
        elif len(m2) is 0:
            match = 'None'
    return match
                    

for infile in listing:
    if str(infile).endswith('-blast.txt'):
        OTU_tax_file = open(infile, 'r')
        print ('Analysing %s' % infile)
        label = str.split(infile, ".txt")[0]
        print label
        outfile = open('%s-new.txt' % label, 'w')
        outfile.write('OTUid\tpath\tSuperkingdom\tKingdom\tSubkingdom\tInfrakingdom\t'
                        'Superphylum\tPhylum\tSubphylum\tInfraphylum\tSuperclass\t'
                        'Class\tSubclass\tInfraclass\tSuperorder\tOrder')        
        match = ''
        OTUid = ''
        path = ''
        path_list = ''
        # Extract taxonomic ranks from MEGAN taxonomic path
        for row in OTU_tax_file: 
            OTUid, path = row.strip().split("\t")
            path_list = path.split(";")
            print path
            outfile.write('\n%s\t%s\t' % (OTUid, path))
            z = len(path_list)-1
            for i in range(z, 0, -1):  # starting from lowest rank in MEGAN path
                t = path_list[i]
                print "path list %i: %s" % (i, t)
                match = get_match(t, tax_ref)
                print "match:"
                if match != 'None':
                    print "match found"
                    break
                else:
                    print "next"
                    next
            if match != 'None':
                #parts = match[0].split('\t')
                #print "parts:"
                #print parts
                #for part in parts:
                #    outfile.write('\t%s' % part)
                outfile.write(match.strip())
            elif match == 'None':
                outfile.write('root\troot\troot\troot\troot\troot\troot\troot\troot\troot\troot\troot\troot\troot')

            if path == 'root':    
                outfile.write('root\troot\troot\troot\troot\troot\troot\troot\troot\troot\troot\troot\troot\troot')

OTU_tax_file.close()                
outfile.close()      
