# The file to initiate and load data 
# depend on ComMA https://github.com/walterxie/ComMA

# names presented in files
input.names <- c("16S","18S","26S","ITS","ShCO1","FolCO1")
# names presented in figures and tables
getOutputNames <- function(input.names) {
  output.names <- gsub("FolCO1", "COI-650", input.names)
  output.names <- gsub("ShCO1", "COI-300", output.names)
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

# get plot names from subplots vector separated by sep
getPlot <- function(subplots, sep="-") {
  sapply(strsplit(as.character(subplots), sep), "[[", 1)
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

###### taxa assignment by reads #####
# "ARCHAEA", "BACTERIA", "CHROMISTA", "PROTOZOA", "FUNGI", "PLANTAE", "ANIMALIA", "EUKARYOTA", "PROKARYOTA", "PROTISTS"
# PROKARYOTA: all prokaryotes (Bacteria + Archaea)
# EUKARYOTA: all eukaryotes
# PROTISTS: "CHROMISTA|PROTOZOA" all micro-eukaryotes
# tt <- getTaxaTable("16S", taxa.group="BACTERIA")
getTaxaTable <- function(data.set=c("16S","18S","26S","ITS","FolCO1","ShCO1"), 
                         data.folder="./data/Taxonomy_tables", taxa.group="all", rank="kingdom") {
  data.set <- match.arg(data.set)
  
  tt.file.path <- file.path(data.folder, paste(data.set, "taxonomy_table.txt", sep="_"))
  require(ComMA)
  taxa.table <- ComMA::readTaxaTable(tt.file.path, matrix.name=data.set, taxa.group=taxa.group, rank=rank)	

  return(taxa.table)
}


###### Trees #####
# tre <- getPhyloTree("16S", "PROKARYOTA")
getPhyloTree <- function(data.set=c("16S","18S","26S","ITS","FolCO1","ShCO1"), 
                         taxa.group="PROKARYOTA", data.folder="./data/Trees", verbose=FALSE) {
  data.set <- match.arg(data.set)
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
    if(verbose) print(env)
  } else {
    env <- NULL
    warning("Cannot find enviornmental data: ", input.f, " !\n") 
  }
  return(env) 
}


######## elevations #######
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




