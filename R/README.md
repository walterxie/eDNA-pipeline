
![Modules](/Modules) folder contains all functions used in the following R code.

## 1. Create Community Matrix and Rarefaction Table

### createAllCommunityMatrix.r
Create community matrix which is the transverse matrix of community matrix in *vegan* package.
Code is tested using output from usearch8.0.1623

* Input: 1) OTUs mapping file "out.up" from USEARCH OTU clustering; 
* Input: 2) OTU representative sequences "otus.fasta" from USEARCH OTU clustering;
* Input: 3) Chimeras sequences chimeras.fasta from Uchime; 
* Input: 4) Duplicate sequences mapping file "derep.uc" from USEARCH de-replication. 
* Input: 5) SraRunTable.txt, a SRA mapping file to map SRA code to subplot name.
* Output: 1) community matrix CSV file, rows are OTUs, columns are samples.  

Note: USEARCH 8 mixies all sample-based data into one pool,
therefore the sample information of each duplicate read would be lost during 
de-replication process. In this code, we retrive the sample of each duplicate 
read from the mapping file derep.uc created by a modified command below to 
create the community matrix retaining the correct reads' distribution of samples.  

```
$USEARCH -derep_fulllength ./qc/denoised.fasta -fastaout ./qc/derep.fasta -sizeout -uc ./qc/derep.uc
```

### createAllDiversitiesOTUsTable.r
Create the rarefaction table and OTU threshold table. This is time-consuming.

* Input: 1) community matrix file, such as 16S.csv; 
* Output: 1) rarefaction table to plot rarefaction curves;  
* Output: 2) OTU threshold table to plot diversities ranging from 90-100% OTU clustering thresholds.  

### createAllRarefactionTable.r
If you do not need clustering through different thresholds, this script is faster to generate the rarefaction table at 97% threshold only.

* Input: 1) community matrix file, such as 16S.csv; 
* Output: 1) rarefaction table to plot rarefaction curves.  

## 2. Community Analysis

### allDiversitiesOTUs.r
Plots of diversities using cutoff thresholds ranging from 90-100% for OTU classification 

* Input: 1) OTU threshold table; 
* Output: 1) the figure of diversities ranging from 90-100% OTU clustering thresholds.  

### allRarefactions.r
Rarefaction curves for diversities estimated.

* Input: 1) rarefaction table; 
* Output: 1) the figure of rarefaction curves.  

### allSampleCount.r
Relative proportion of OTUs inferred by read count.

* Input: 1) community matrix file, such as 16S.csv; 
* Output: 1) the figure of relative proportion of OTUs inferred by read count;
* Output: 2) the figure of OTUs across the sites.

### allStatistics.r
Generate statistical summary of communities regarding sequences and OTUs.

* Input: 1) community matrix file, such as 16S.csv; 
* Output: 1) Latex tables, statistical summary regarding sequences and OTUs.

### allTaxonomyPhylum.r
The percentage of OTUs assigned to phyla. 

* Input: 1) a taxonomic identification file from MEGAN but modified by expert, such as taxonomy97phyla.txt; 
* Output: 1) the figure of OTUs assigned to phyla.

### allWithinBetweenPlots.r
Box and whisker plots of turnover (normalised pairwise effective beta diversity) within (red) and between plots (blue) 

* Input: 1) community matrix file, such as 16S.csv; 
* Output: 1) the figure of turnover within and between plots.

### allMDSBySubplots.r
Multivariate ordination of samples using non-metric multidimensional scaling of effective beta diversity for paired subplots.

* Input: 1) community matrix file, such as 16S.csv;
* Output: 1) the figure of MDS for paired subplots. 

### allGeneCorrolation.r
Pairwise community matrix correlations of effective beta diversity between data sets.

* Input: 1) community matrix file, such as 16S.csv; 
* Output: 1) the figure of MDS for correlations;
* Output: 2) the table of correlations and significances.

### allMaxDiv.r
The probability of having maximum gamma diversity or maximum effective beta diversity of all possible combinations of 4 plots. 
The number of plots (4) can be changed in the code.

* Input: 1) community matrix file, such as 16S.csv; 
* Output: 1) the table of probabilities.

### allMaxDivHeatmap.r
Create heat-map for *allMaxDiv.r*.

* Input: 1) the table of probabilities; 
* Output: 1) the heat-map.

### allMaxRemainedDiversity.r
Maximum remained gamma and effective beta diversity as a function of number of sites. 
1 is the most important and removed at the last, 10 is the least important and removed in the beginning.

* Input: 1) community matrix file, such as 16S.csv; 
* Output: 1) the table of ranks;
* Output: 2) the table of Spearman correlations and their significances.

### avgRanks.r
Extension code of *allMaxRemainedDiversity.r*.

* Input: 1) the table of ranks; 
* Output: 1) the table of means and standard deviations of given ranks.


## 3. Take Environmental Variables Into Account 

### allElevationAlpha.r
Regression of effective alpha diversity as a function of elevation.

* Input: 1) community matrix file, such as 16S.csv; 
* Output: 1) the figure of regression lines;
* Output: 2) the table of Mantel's test result.

### allElevationDiversitiesByPlots.r
Regression of effective beta diversity and difference in elevation.

* Input: 1) community matrix file, such as 16S.csv; 
* Output: 1) the figure of regression lines;
* Output: 2) the table of Mantel's test result.

### allRedundancyAnalysis.r
Constrained ordination of community data with environmental data as constraining variables, was carried out using the *[capscale](http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/capscale.html)* function, 
a non-Euclidean generalization of redundancy analysis, from the R package *vegan*. 
Three ordination scenarios were tested: 

  i) Models were constructed containing each of the fifteen environmental variables in isolation. 
  
  ii) A model was constructed using the combined set of variables with VIF < 10. 
  
  iii) More conservative models were constructed by using subsets of the variables with VIF < 10 chosen 
  by stepwise forward and backward selection model building procedures.

* Input: 1) community matrix file, such as 16S.csv;
* Input: 2) table of environmental location and conditions for sampling sites, such as plot_elevations.txt;
* Input: 3) soil chemistry, such as LJ12027.txt; 
* Output: 1) the figure of correlation plots of environmental variables;
* Output: 2) the figure of models in ii); 
* Output: 3) the figure of models using forward selection in iii);
* Output: 4) the figure of models using backward selection in iii);
* Output: 5) the result table of distance-based redundancy analysis and their ANOVA tests in each step;
* Output: 6) the table of constrained and unconstrained inertia changes during analysis.
