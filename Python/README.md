# eDNA-pipeline

In development

## Python Scripts

1. phylogenetic\_alpha\_v3.py is to calculate phylogenetic diversity in various forms. 
It is a bit slow with big trees of miseq data (especially calculating Rao's "Q"). 
  
2. match\_new\_vs\_old_taxonomy.py is used to convert old NCBI taxonomy to 
the new taxonomy, and make taxonomy tables from taxa paths, such as MEGAN-exported paths.
Note: The last ";" in the taxa path has to be removed. 
  
3. New\_taxonomy\_from\_PLOSONE is the new taxonomy scheme that comes from:
Ruggerio et al. 2015, A Higher Level Classification of All Living Organisms. 
http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0119248

4. make\_taxonomy\_table\_v5.py