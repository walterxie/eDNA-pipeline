# -*- coding: utf-8 -*-
"""
Makes a table of OTU counts from UPARSE algorithm output and dereplication .uc file.

@author: A Dopheide
"""

def make_otutable(label, up_file, derep_file):
    print("%s: Processing uparse file..." % label)
    samples = []
    otus = []
    otu_matches = {}
    chimeras = {}

    for row in open(up_file, 'r'):
        fields = row.strip().split("\t")
        query = fields[0]
        query_id = query.split(";")[0]
        #print(query_id)
        sample = query_id.split("|")[2]
        if sample not in samples:
            samples.append(sample)
        best_match = fields[4]
        if "otu" in fields[1]:  # Query is OTU centroid
           # print("OTU centroid")
            otus.append(query_id)
            otu_matches[query_id] = query_id
        elif "match" in fields[1]:  # Query is 97 % OTU match
            #print("match")
            otu_matches[query_id] = best_match
        elif "chimera" in fields[1]:  # Query is a chimera
            #print("chimera")
            chimeras[query_id] = best_match

    # Set up bins and tables for otu counts
    bins = [[0 for row in range(len(samples))] for col in range(len(otus))]
    otutable = open(("%s_otutable.txt" % label), "w")

    # Save chimeras to file
    chimeras_out = open(("%s_chimeras.txt" % label), "w")
    print("Writing chimeras...")
    for chimera, best_match in chimeras.iteritems():
        chimeras_out.write("%s\t%s\n" % (chimera, best_match))
    chimeras_out.close()

    # Write otu counts from uparse out file to otu table
    for query_id, match in otu_matches.iteritems():
        sample = query_id.split("|")[2]
        bins[otus.index(match)][samples.index(sample)] += 1

    # Get replicate sequence counts from derep.uc file, add to bins
    print("%s: Processing derep.uc file...")
    for row in open(derep_file, 'r'):
        fields = row.strip().split("\t")
        query = fields[8]
        sample = query.split("|")[2]
        if fields[0] is "S":  # Unique sequence/seed
            #print(query)
            #print("derep seed")
            otu = query
        elif fields[0] is "H":  # Replicate sequence/hit
            #print(query)
            #print("derep hit")
            derep_hit = fields[9]  # Get matching seed for hit
            otu = otu_matches.get(derep_hit, "Not found")  # Get matching otu centroid for seed
            if otu is not "Not found":  # i.e. isn't a chimera
                bins[otus.index(otu)][samples.index(sample)] += 1

    for item in samples:
        otutable.write("\t%s" % item)
    n = 0
    for row in bins:
        otutable.write("\n%s" % otus[n])
        for item in row:
            otutable.write("\t%s" % item)
        n += 1

    otutable.close()
    up_file.close()
    derep_file.close()
    print("Finished %s" % label)
