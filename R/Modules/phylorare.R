## R function for the calculation of rarefaction of phylogenetic diversity of multiple sites.
## by David Nipperess, Macquarie University, Australia (david.nipperess@mq.edu.au)
## Last edited 21st January 2014.

## This software is made available free of charge.
## No technical support for using this software will be provided.
## You must download and install the ape package (available from http://cran.r-project.org/) to use this function.

## ---------------

phylorare <- function (x, phy, m=1, subsampling = "individual", replace = F) {

# where x is a community data table (as in the vegan package) with species/OTUs as columns and sites as rows. Columns are labelled with the names of the species/OTUs. Rows are labelled with the names of the sites. Data is either abundance or incidence.
# where phy is a phylogenetic tree stored as a phylo object (as in the ape package) with terminal nodes labelled with names matching those of the community data table. Note that the function trims away any terminal taxa not present in the community data table, so it is not necessary to do this beforehand.
# where m is the number of individuals, sites or species to which x is to be rarefied.
# where subsampling indicates whether the subsampling will be by "individual" (default), "site" or "species". When there are multiple sites, rarefaction by individuals or species is done for each site separately.
# where replace = T indicates that probabilities are treated as if populations of individuals, sites or species are of infinite size.

require(ape)

### step 1: trimming the tree to match the community data table thus creating a "community tree" (sensu Cam Webb).

if (length(phy$tip.label) > length(x[1,])) {
	phy <- drop.tip (phy, which(!phy$tip.label %in% colnames(x))) }

# script is modified from that of Paradis (2006) "Analysis of Phylogenetics and Evolution with R", Springer.

### step 2: converting a community tree into a MRP matrix

# A MRP matrix, used in supertree methods, is where the membership of an OTU in a clade spanned by a branch is indicated by a 0 (no) or 1 (yes). Unlike supertree MRP matrices, our matrix includes terminal branches.
# the new phylo object model for ape 2.0 broke the original code. The following replaces that code.

phylomatrix <- matrix (0, length(phy$tip.label), length(phy$edge.length))
for (i in 1:length(phy$tip.label)) {
	lineage <- which (phy$edge[,2] == i)
	node <- phy$edge[lineage,1]
	while (node > length(phy$tip.label)+1) {
		branch <- which (phy$edge[,2] == node)
		lineage <- c (lineage, branch)
		node <- phy$edge[branch,1]
		}
	phylomatrix[i,lineage] = 1
	}

# this script fills a matrix of 0's with 1's corresponding to the incidence of an OTU on a particular branch.
# the code is pretty slow on large trees. I will try timing it for the mammal supertree.

### step 3: re-ordering the OTUs of the occurrence table and MRP matrix to match.

phylomatrix <- phylomatrix[sort.list(phy$tip.label), ]
x <- x[ ,sort.list(colnames(x))]

# data are sorted to a common ordering standard, that is alphabetic order, so that OTUs match up.

### step 4: creating a community phylogeny matrix from a MRP matrix and an occurrence matrix

x <- as.matrix (x)
if(subsampling=="species") {
	x <- ifelse(x>0,1,0)}
	
commphylo <- x %*% phylomatrix  # commphylo will now be either the no. of individuals or species per branch per sample

if(subsampling=="site") {
	commphylo <- ifelse (commphylo > 0, 1, 0) # convert to incidence form
	commphylo <- colSums(commphylo) # reduce to vector of branch frequencies
	}

# the above code performs simple matrix multiplication to produce a composite branch by sample matrix (the community phylogeny matrix) where each cell now contains total abundance for each branch in each sample.

### step 5: convert above matrix into probabilities based on size of m.

if(subsampling=="site") {
	N <- length(x[,1]) # this gives the no. of sites
	}
else {
	N <- rowSums(x) # this gives either the no. of individuals or the no. of species
	}

if (replace == T) {
		p <- 1-(1-(commphylo/N))^m
		}
else {
		p <- 1-(exp(lchoose((N-commphylo),m)-lchoose(N,m)))
		}

### step 6: calculating rarefied PD per sample

pdrare <- p %*% phy$edge.length

return (pdrare)

}