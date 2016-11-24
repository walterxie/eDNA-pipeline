# -*- coding: utf-8 -*-
"""
This takes a file of taxonomic paths exported from MEGAN metagenome 
analyzer, matches up the taxonomic ranks above order level with taxonomy from 
Ruggerio et al. 2015, and below order level with taxonomy from NCBI, and 
outputs a table of OTU ids, taxconomic paths, and each rank down to genus or
whichever rank is the lowest level (usually a higher level than genus).

@author: A Dopheide
"""
import os
#import re
import glob

rank1 = set()
rank2 = set()
domains = set()
majorclades = set()
superkingdoms = set()
kingdoms = set()
subkingdoms = set()
superphyla = set()
phyla = set()
subphyla = set()
superclasses = set()
classes = set()
subclasses = set()
infraclasses = set()
superorders = set()
orders = set()
suborders = set()
superfamilies = set()
families  = set()
subfamilies = set()
genera = set()
everything = set()

### Compile sets of taxa ranks
### Use NCBI taxonomy for orders and below
### First get orders through subfamilies
os.chdir('../Taxonomy/')
os.chdir('G:/Documents\PhD/taxonomy_test')
with open('Database_stuff/NCBI_taxa_list.txt', 'r') as tax_ref_NCBI:
    for row in tax_ref_NCBI:
        fields = row.strip().split('\t')
        node = fields[1]    
        rank = fields[3]
        #print(i)        
        print(node)
        print(rank)            
        if rank != "species" and rank != "genus":
            if node not in everything:
                everything.add(node)
        if rank == "order":
            orders.add(node)
        elif rank == "suborder":
            suborders.add(node)
        elif rank == "superfamily":
            superfamilies.add(node)
        elif rank == "family":
            families.add(node)
        elif rank == "subfamily":
            subfamilies.add(node)

### Now get genera
with open('Database_stuff/NCBI_taxa_list.txt', 'r') as tax_ref_NCBI:     
    for row in tax_ref_NCBI:
        nodebits = list()
        fields = row.strip().split('\t')
        node = fields[1]    
        rank = fields[3]
        if rank == "genus":
            print(node) 
            if node not in everything: # Ensures genus isn't present among higher-level ranks
                genera.add(node)
        elif rank == "species": # Some species aren't represented among genera
            print(node) 
            g = node.split(' ')[0] # Gets genus from species name 
            if g not in everything: # If genus not already present
                genera.add(g) # Add to pool of genus names

### Remove some dubious genus names, e.g. "Actinobacteridae"
len(genera) # 77753 genera
genera = [g for g in genera if not g.endswith('idae') and not g.endswith('aceae')]
len(genera) # now 77708 genera
          
### Get high-level taxonomic rank paths from Ruggiero et al. 2015 taxonomy
### (Note this omits a few groups, such as Cryptomycota)
### Compile set of taxonomic paths and lowest rank names 
tax_paths = {}
l = []
n = 1
with open('Database_stuff/New_taxonomy_from_PLOSONE_2015_cleaned_up.txt', 'r') as tax_ref_high_level:
    for row in tax_ref_high_level:
        print (n)
        if n is 1:
            n += 1 # Skip the first row (rank names)
        else:
            fields = row.strip('\n').split('\t')
            last = fields[13] # Order or equivalent lowest rank level
            print (last)
            tax_paths[last] = row
            l.append(last)
            n += 1

### Function to find matching taxonomic path for each taxon
def get_match(t, tax_paths):
    m = [paths for taxon, paths in tax_paths.items() if t.lower() == taxon.lower()]
    #print (len(m))
    if len(m) is 1:  # exact match
        match = m[0]
    elif len(m) > 1: # multiple matches (e.g. Chlorophyta and CHLOROPHYTA)
        match = m[0] # take the first match
    elif len(m) is 0: # no exact match
        ### So check for less exact matches? 
        ### Accounts for some name differences, but also results in some incorrect ones
        #m2 = [paths for taxon, paths in tax_paths.items() if t.lower() in taxon.lower()]            
        #if len(m2) is 1:
            #match = m2[0]
        #elif len(m2) > 1:
            #match = m2[0]
        #elif len(m2) is 0: 
        match = 'None' # no match (so try next rank up)
    return match

### Covert taxonomic paths to table of ranks for each gene
os.chdir('../Data/')
genelist = ['16S','18S','26S','ITS','FolCO1','ShCOI']
for g in genelist:
    files = glob.glob('%s_*OTUs_nt_paths.txt' % g) # Taxonomic paths file
    with open(files[0], 'r') as OTU_paths_file:
        with open('%s_taxonomy_table.txt' % (g), 'w') as outfile:
            outfile.write('OTUid\tpath\tSuperkingdom\tKingdom\tSubkingdom\tInfrakingdom\t'
                                'Superphylum\tPhylum\tSubphylum\tInfraphylum\tSuperclass\t'
                                'Class\tSubclass\tInfraclass\tSuperorder\tOrder\tSuborder\t'
                                'Superfamily\tFamily\tSubfamily\tGenus') 
            
            for row in OTU_paths_file:
                match = ''
                order = ''
                suborder = ''
                superfamily = ''
                family = ''
                subfamily = ''
                genus = ''
                OTUid, path = row.strip().split('\t')
                path_list = path.strip(';').split(';') # Split path into ranks
                outfile.write('\n%s\t%s' % (OTUid, path))

                ### Get Order and higher ranks from Ruggiero et al. 2015 taxonomy
                #print(OTUid)
                #print(path)
                z = len(path_list)
                for i in range(z, 0, -1): # Starting from lowest rank in MEGAN path...
                    t = path_list[i-1] # Get the last (or next lowest) rank
                    #print ("path list %i: %s" % (i-1, t))
                    match = get_match(t, tax_paths).strip() # Get matching high-level ranks
                    if match is not 'None': # Matching high-level ranks path found
                        print ("matching path: %s" % match)
                        #print ("match found")
                        break 
                    else: # If match not found, move on to next lowest rank
                        next 
                outfile.write('\t%s' % match) # Write high-level ranks to file
                ### Get lowest rank in matching high-level taxonomy, 
                ### carry over to lower ranks if other lower ranks not found
                parts = match.strip().split('\t')
                last = parts[len(parts)-1] 
                ### Organise orders and below according to NCBI taxonomy        
                ### i.e. check if the lowest rank name in the path is an order,
                ### suborder, family (etc.), and organise accordingly    
                order = next((d for d in path_list if d in orders), last)
                suborder = next((d for d in path_list if d in suborders), order)
                family = next((d for d in path_list if d in families), suborder)
                subfamily = next((d for d in path_list if d in subfamilies), family)
                genus = next((d for d in path_list if d in genera), subfamily)            
                print('order: %s' % order)
                print('suborder: %s' % suborder)
                print('family: %s' % family)
                print('subfamily: %s' % subfamily)
                print('genus: %s' % genus)
                outfile.write(("\t%s\t%s\t%s\t%s\t%s") % (order, suborder, family, subfamily, genus))
