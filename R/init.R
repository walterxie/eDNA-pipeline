# The file to initiate and load data 
# depend on ComMA https://github.com/walterxie/ComMA

# names presented in files
input.names <- c("16S","18S","26S","ITS","ShCO1","FolCO1")
# names presented in figures and tables
getOutputNames <- function(input.names) {
  output.names <- gsub("FolCO1", "COI-650", input.names, ignore.case = T)
  output.names <- gsub("ShCO1", "COI-300", output.names, ignore.case = T)
}
# convert taxa names into valid name in taxa table file
getTaxaNames <- function(taxa.names) {
  taxa.names <- gsub("protists", "CHROMISTA|PROTOZOA", taxa.names, ignore.case = T)
  taxa.names <- gsub("animals", "ANIMALIA", taxa.names, ignore.case = T)
  taxa.names <- gsub("all OTUs", "all", taxa.names, ignore.case = T)
}
# find the rank associated to taxa group to be analysed
getRank <- function(taxon) {
  if (toupper(taxon) == "PROKARYOTA" || toupper(taxon) == "EUKARYOTA") "superkingdom"
  else "kingdom"
}

# most abundant 150 OTUs
most.abundant.OTU = 150

# add postfix for various names
postfix <- function(..., sep, min2=TRUE, by.plot=FALSE) {
  name <- paste(..., sep = sep)
  if (min2) 
    name <- paste(name, "min2", sep = sep)
  if (by.plot) 
    name <- paste(name, "by", "plot", sep = sep) 
  return(name)
}

######## load community matrix #######
# by_plot=T plot based, min2=T no singleton 
# cm <- getCommunityMatrix("16S")
getCommunityMatrix <- function(data.set=c("16S","18S","26S","ITS","FolCO1","ShCO1","Vegetation"), 
                               min2=TRUE, by.plot=FALSE, data.folder="./data/OTU_tables") {
  data.set <- match.arg(data.set)

  if (data.set=="Vegetation") {
    if (!by.plot)
      stop("Vegetation only has plot based community matrix !")
    cm.file.path <- file.path(data.folder, "..", "LBI_Trees_Saplings_SBA.csv")
  } else {
    fn <- postfix(data.set, "otutable", sep="_", min2=min2, by.plot=by.plot)
    # e.g. data/16S_otutable.txt
    cm.file.path <- file.path(data.folder, paste(fn, "txt", sep="."))
  }
  print(cm.file.path)
  require(ComMA)
  # always set minAbund=1 here
  community.matrix <- ComMA::readCommunityMatrix(cm.file.path, matrix.name=data.set, minAbund=1)
  
  return(community.matrix)
}
# t.cm <- transposeCM(cm)

# taxa.table <- getTaxaTable("16S", taxa.group="assigned")
# tt.sub <- subsetTaxaTable(taxa.table, taxa.group="BACTERIA", rank="kingdom")
# cm.taxa <- mergeCMTaxa(cm, taxa.table)
# taxa.assign <- assignTaxaByRank(cm.taxa)
# sum(taxa.assign[[rank]])

# cm.taxa <- ComMA::subsetCM(cm, tt, taxa.group="BACTERIA", rank="kingdom")

getIdentifiedCM <- function(data.set=c("16S","18S","26S","ITS","FolCO1","ShCO1"), 
                            by.plot=FALSE, min2=TRUE, 
                            cm.folder="./data/OTU_tables", tt.folder="./data/Taxonomy_tables") {
  data.set <- match.arg(data.set)
  if (data.set == "16S")
    input.tg <- "PROKARYOTA"
  else 
    input.tg <- "EUKARYOTA"
  # no singleton
  cm <- getCommunityMatrix(data.set, min2=min2, by.plot=by.plot, data.folder=cm.folder)
  tt <- getTaxaTable(data.set, taxa.group=input.tg, rank="superkingdom", data.folder=tt.folder)
  # ITS tt only up to "order"
  return(ComMA::subsetCM(cm, tt, col.ranks=c("kingdom", "phylum", "class", "order")))
}

###### taxa assignment by reads #####
# "ARCHAEA", "BACTERIA", "CHROMISTA", "PROTOZOA", "FUNGI", "PLANTAE", "ANIMALIA", "EUKARYOTA", "PROKARYOTA", "PROTISTS"
# PROKARYOTA: all prokaryotes (Bacteria + Archaea)
# EUKARYOTA: all eukaryotes
# PROTISTS: "CHROMISTA|PROTOZOA" all micro-eukaryotes
# tt <- getTaxaTable("16S", taxa.group="BACTERIA")
getTaxaTable <- function(data.set=c("16S","18S","26S","ITS","FolCO1","ShCO1"), 
                         taxa.group="all", rank="kingdom", data.folder="./data/Taxonomy_tables") {
  data.set <- match.arg(data.set)
  
  tt.file.path <- file.path(data.folder, paste(data.set, "taxonomy_table.txt", sep="_"))
  require(ComMA)
  taxa.table <- ComMA::readTaxaTable(tt.file.path, matrix.name=data.set, taxa.group=taxa.group, rank=rank)	

  return(taxa.table)
}

getTaxaRef <- function(data.folder="./data") {
  taxa.ref <- ComMA::readFile(file.path(data.folder, "New_taxonomy_from_PLOSONE_2015_cleaned_up.txt"), 
                             row.names = NULL, quote = "", comment.char = "")
  # use samll case in ComMA
  colnames(taxa.ref) <- tolower(colnames(taxa.ref))
  return(taxa.ref)
}

# programmatically get sub-dataset 
getCommunityList <- function(genes=c("16S","18S","26S","ITS","FolCO1","ShCO1"),
                             genes.taxa=list(list("16S","bacteria"),list("18S","protists"),list("18S","fungi"),
                                             list("18S","animals"),list("26S","fungi"),list("ShCO1","animals")), 
                             by.plot=TRUE, col.ranks=c("superkingdom", "kingdom"), drop.taxa=TRUE ) {
  taxa=c("all", "assigned", "ARCHAEA", "BACTERIA", "CHROMISTA", "PROTOZOA", 
         "CHROMISTA|PROTOZOA", "FUNGI", "PLANTAE", "ANIMALIA", "EUKARYOTA")
  
  # data frame for statistics
  cm.taxa.list <- list()
  require(ComMA)
  for (z in genes.taxa) {
    # no singletons
    min2=TRUE
    gene <- getOutputNames(z[[1]])
    taxon <- getTaxaNames(z[[2]])
    rank <- getRank(taxon)
    if ( !tolower(z[[1]]) %in% tolower(genes) || !tolower(taxon) %in% tolower(taxa) )
      stop("Invalid name in genes.taxa ! ", gene, " (", z[[1]], ") or ",  taxon, " (", z[[2]], ")")
    
    cat("\n", gene, "data set", taxon, "group at", rank, ",", ifelse(min2, "exclude", "include"), 
        "singletons, samples are based on", ifelse(by.plot, "plot", "subplot"), ".\n") 
    
    cm <- getCommunityMatrix(z[[1]], min2=min2, by.plot=by.plot)
    tt <- getTaxaTable(z[[1]], taxa.group=taxon, rank=rank)
    
    cm.taxa <- ComMA::mergeCMTaxa(cm, tt, has.total = 0, col.ranks = col.ranks, preprocess = F)
    
    if (drop.taxa)
      cm.taxa <- cm.taxa[,-which(names(cm.taxa) %in% col.ranks)]
    
    cm.taxa.list[[paste(gene, z[[2]])]] <- cm.taxa
  }
  cat("\n")
  
  return(cm.taxa.list)
}


###### Trees #####
# tre <- getPhyloTree("16S", "PROKARYOTA")
getPhyloTree <- function(data.set=c("16S","18S","26S","ITS","FolCO1","ShCO1"), 
                         taxa.group=NA, data.folder="./data/Trees", verbose=FALSE) {
  data.set <- match.arg(data.set)
  # set taxa.group=NA to auto choose taxa.group
  if (is.na(taxa.group)) {
    if (data.set == "16S")
      taxa.group <- "PROKARYOTA"
    else 
      taxa.group <- "EUKARYOTA"
  }
  
  input.f <- file.path(data.folder, paste(data.set, tolower(taxa.group), "min2.tre", sep = "-"))
  if (file.exists(input.f)) {
    require(ape)
    cat("\nLoad tree from", input.f, "\n") 
    tree <- read.tree(input.f)
    if(verbose) print(tree)
  } else {
    tree <- NULL
    warning("Cannot find tree file: ", input.f, " !\n") 
  }
  return(tree)
}

######## meta data of samples #######
# env.plot <- getEnvData
getEnvData <- function(by.plot=TRUE, data.folder="./data/Environmental_data", verbose=FALSE) {
  input.f <- file.path(data.folder, paste("LBI_all_env_data_by", 
                       ifelse(by.plot, "plot.txt", "subplot.txt"), sep = "_"))
  if (file.exists(input.f)) {
    require(ComMA)
    cat("\nLoad enviornmental data from", input.f, "\n") 
    env <- ComMA::readFile(input.f)
    #if (by.plot) env$Plot <- rownames(env)
    if(verbose) print(env)
  } else {
    env <- NULL
    warning("Cannot find enviornmental data: ", input.f, " !\n") 
  }
  env[is.na(env)] <- "Unknown"
  return(env) 
}

### Load between plot distances ###
getBetweenPlotDistance <- function(data.folder="./data/Environmental_data", verbose=TRUE) { 
  input.f <- file.path(data.folder, "LBI_between-plot_distances.txt")
  if (file.exists(input.f)) {
    require(ComMA)
    if(verbose)
      cat("\nLoad between plot distances from", input.f, "\n") 
    plot.dists <- ComMA::readFile(input.f, verbose=verbose)
    plot.dists$pair <- paste(plot.dists$plot1, plot.dists$plot2)
    plot.dists$pair <- plotSort(plot.dists$pair)
  } else {
    plot.dists <- NULL
    warning("Cannot find between plot distances: ", input.f, " !\n") 
  }
  return(plot.dists) 
}

### Load within plot distances ###
getWithinPlotDistance <- function(data.folder="./data/Environmental_data", verbose=TRUE) { 
  input.f <- file.path(data.folder, "LBI_subplot_distances.txt")
  if (file.exists(input.f)) {
    require(ComMA)
    if(verbose)
      cat("\nLoad within plot distances from", input.f, "\n") 
    s <- ComMA::readFile(input.f, verbose=verbose)
    subplot.dists <- data.frame(t(combn(names(s), 2)), dist=s[lower.tri(s)])
    subplot.dists$pair <- paste0(subplot.dists$X1, subplot.dists$X2)
    subplot.dists$pair <- subplotSort(subplot.dists$pair)
  } else {
    subplot.dists <- NULL
    warning("Cannot find between plot distances: ", input.f, " !\n") 
  }
  return(subplot.dists) 
}


######## plot vs subplots #######
# get plot names from subplots vector separated by sep
getPlot <- function(subplots, sep="-") 
  sapply(strsplit(as.character(subplots), sep), "[[", 1)

### sort pairs of plot names ###
plotSort <- function(x)
  sapply(lapply(strsplit(x, " "), sort), paste, collapse=" ")
### Sort subplot pairs ###
subplotSort <- function(x)
  sapply(lapply(strsplit(x, NULL), sort), paste, collapse="")







######## elevations dist #######
getElevPlotDist <- function(plot.names, env.byplot) { 
  colElev = 1
  # case insensitive
  matched.id <- match(tolower(plot.names), tolower(rownames(env.byplot)))
  matched.id <- matched.id[!is.na(matched.id)]
  # match 
  env.plot.match <- env.plot[matched.id, ]
  
  cat("Find", nrow(env.plot.match), "plots having elevations, community matrix has", 
      length(plot.names), "plots, meta-data file has", nrow(env.plot), "plots.\n")
  
  return(dist(env.plot.match[,colElev]))
}

###### table to plot Phylo Rarefaction ##### 
getPhyloRareTable <- function(expId, isPlot, min2, taxa.group="assigned") {
  n <- length(matrixNames) 
  # hard code for Vegetation that only has plot and always keep singletons
  if (expId == n) {
    mid.name <- postfix("all", TRUE, FALSE, sep="-")
  } else {
    mid.name <- postfix(taxa.group, isPlot, min2, sep="-") 
  }
  
  input.f <- file.path(workingPath, "data", "pdrf", paste(matrixNames[expId], mid.name, "phylorare", "table.csv", sep="-"))
  if (file.exists(input.f)) {
    phylo.rare.df <- read.csv(file=input.f, head=TRUE, sep=",", row.names=1, check.names=FALSE)
    if(verbose) 
      cat("\nUpload phylo rarefaction table from", input.f, "\n") 
  } else {
    phylo.rare.df <- NULL
    cat("Warning: cannot find phylo rarefaction table", input.f, "\n") 
  }
  phylo.rare.df
}

###### table to plot Rarefaction ##### 
getRarefactionTableTaxa <- function(expId, isPlot, min2, taxa.group, div="alpha1") {
  pathFileStem <- file.path(workingPath, "data", "rf", paste(matrixNames[expId], 
                    postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-"))
  input.f <- paste(pathFileStem, "rare", div, "table.csv", sep = "-")
  if (file.exists(input.f)) {
    rare.df <- read.csv(file=input.f, head=TRUE, sep=",", row.names=1, check.names=FALSE)
    if(verbose) 
      cat("\nUpload rarefaction table per sample from", input.f, "\n") 
  } else {
    rare.df <- NULL
    cat("Warning: cannot find rarefaction table per sample", input.f, "\n") 
  }
  rare.df
}

getRarefactionTable <- function(expId, isPlot, min2) {
  n <- length(matrixNames) 
  matrixName <- matrixNames[expId]
  # hard code for Vegetation that only has plot and always keep singletons
  if (expId == n) {
    matrixName <- postfix(matrixName, TRUE, FALSE, sep="-")
  } else {
    matrixName <- postfix(matrixName, isPlot, min2, sep="-") 
  }
  
  inputRDT <- file.path(workingPath, "data", paste(matrixName, "rarefaction-table.csv", sep="-"))
  if(verbose) 
    cat("\nUpload rarefaction table : from", inputRDT, "\n") 
  
  rarefactionTable <- read.csv(file=inputRDT, head=TRUE, sep=",", row.names=paste(levels, qs, sep=""), check.names=FALSE)
}

###### dissimilarity matrix #####
# Dissimilarity matrix of paired samples
# diss.fun = "beta1-1", "jaccard", "horn.morisita"
getDissimilarityMatrix <- function(expId, isPlot, min2, diss.fun="beta1-1", taxa.group="all") {
  n <- length(matrixNames) 
  # hard code for Vegetation that only has plot and always keep singletons
  if (expId == n) {
    fname <- paste(matrixNames[expId], postfix("all", TRUE, FALSE, sep="-"), diss.fun, sep = "-")
  } else {
    fname <- paste(matrixNames[expId], postfix(taxa.group, isPlot, min2, sep="-"), diss.fun, sep = "-") 
  }
  
  inputB <- file.path(workingPath, "data", "dist", paste(fname, "csv", sep = "."))
  if(verbose) 
    cat("\nUpload", diss.fun, "matrix of", taxa.group, "taxa group(s) from", inputB, "\n") 
  
  diss.matrix <- readFile(file=inputB, sep=",")
  
  return(diss.matrix)
}

###### table to max remained diversity ##### 
getMaxRemainedDiversity <- function(lev.q, taxa.group="assigned") {
  input.f <- file.path(workingPath, "data", "maxrd", paste("max-div", lev.q, taxa.group,"table.csv", sep = "-"))
  if (file.exists(input.f)) {
    max.rd <- read.csv(file=input.f, head=TRUE, sep=",", row.names=1, check.names=FALSE)
    if(verbose) 
      cat("\nUpload max remained diversity table from", input.f, "\n") 
  } else {
    max.rd <- NULL
    cat("Warning: cannot find max remained diversity table", input.f, "\n") 
  }
  max.rd
}

getMaxRemainedDiversityRank <- function(lev.q, taxa.group="assigned") {
  input.f <- file.path(workingPath, "data", "maxrd", paste("max-div-rank", lev.q, taxa.group,"table.csv", sep = "-"))
  if (file.exists(input.f)) {
    max.rd <- read.csv(file=input.f, head=TRUE, sep=",", row.names=1, check.names=FALSE)
    if(verbose) 
      cat("\nUpload max remained diversity rank table from", input.f, "\n") 
  } else {
    max.rd <- NULL
    cat("Warning: cannot find max remained diversity rank table", input.f, "\n") 
  }
  max.rd
}




