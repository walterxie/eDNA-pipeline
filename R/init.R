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
  taxa.names <- gsub("plants", "PLANTAE", taxa.names, ignore.case = T)
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

# programmatically get subsets, drop.taxa TRUE to only return CM,
# genes and taxa.group are used for validation, genes.taxa to choose the subsets. 
getCommunityList <- function(genes=c("16S","18S","26S","ITS","FolCO1","ShCO1"),
                             taxa.group=c("all","assigned","ARCHAEA","BACTERIA","CHROMISTA","PROTOZOA",  
                                    "CHROMISTA|PROTOZOA","FUNGI","PLANTAE","ANIMALIA","EUKARYOTA","PROKARYOTA"),
                             genes.taxa=list(list("16S","prokaryota"),list("18S","eukaryota"),list("26S","eukaryota"),
                                             list("ITS","eukaryota"),list("ShCO1","eukaryota"),list("FolCO1","eukaryota")), 
                             by.plot=TRUE, col.ranks=c("superkingdom", "kingdom"), drop.taxa=TRUE, 
                             pre.cm=TRUE, rm.samples=c(), min.abund=5, mean.abund.thr=0.025 ) {
  # data frame for statistics
  cm.taxa.list <- list()
  require(ComMA)
  for (z in genes.taxa) {
    # no singletons
    min2=TRUE
    gene <- getOutputNames(z[[1]])
    taxon <- getTaxaNames(z[[2]])
    rank <- getRank(taxon)
    if ( !tolower(z[[1]]) %in% tolower(genes) || !tolower(taxon) %in% tolower(taxa.group) )
      stop("Invalid name in genes.taxa ! ", gene, " (", z[[1]], ") or ",  taxon, " (", z[[2]], ")")
    
    cat("\n", gene, "data set", taxon, "group at", rank, ",", ifelse(min2, "exclude", "include"), 
        "singletons, samples are based on", ifelse(by.plot, "plot", "subplot"), ".\n") 
    
    cm <- getCommunityMatrix(z[[1]], min2=min2, by.plot=by.plot)
    
    if (pre.cm)
      cm <- ComMA::preprocessCM(cm, rm.samples=rm.samples, min.abund=min.abund, mean.abund.thr=mean.abund.thr)
    
    if (taxon == "all") {
      # if all, then no merge
      cm.taxa.list[[paste(gene, z[[2]])]] <- cm
    } else {
      tt <- getTaxaTable(z[[1]], taxa.group=taxon, rank=rank)
      cm.taxa <- ComMA::mergeCMTaxa(cm, tt, has.total = 0, col.ranks = col.ranks, preprocess = F)
      
      if (drop.taxa)
        cm.taxa <- cm.taxa[,-which(names(cm.taxa) %in% col.ranks)]
      
      cm.taxa.list[[paste(gene, z[[2]])]] <- cm.taxa
    }
  }
  cat("\n")
  return(cm.taxa.list)
}

# get list of trees, same to getCommunityList
getTreeList <- function(genes=c("16S","18S","26S","ITS","FolCO1","ShCO1"),
                        taxa.group=c("all","assigned","ARCHAEA","BACTERIA","CHROMISTA","PROTOZOA",  
                                     "CHROMISTA|PROTOZOA","FUNGI","PLANTAE","ANIMALIA","EUKARYOTA","PROKARYOTA"),
                        genes.taxa=list(list("16S","prokaryota"),list("18S","eukaryota"),list("26S","eukaryota"),
                                        list("ITS","eukaryota"),list("ShCO1","eukaryota"),list("FolCO1","eukaryota")) ) {
  # data frame for statistics
  tre.list <- list()
  for (z in genes.taxa) {
    gene <- getOutputNames(z[[1]])
    tre <- getPhyloTree(z[[1]], z[[2]])
    tre.list[[paste(gene, z[[2]])]] <- tre
  }
  cat("\n")
  return(tre.list)
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
  if (taxa.group=="all") {
    taxa.group <- "not-assigned"
  } else if (taxa.group=="animals") {
    taxa.group <- "animalia"
  } else if (taxa.group=="plants") {
    taxa.group <- "plantae"
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

######## save figures #######
# use varible name to determine which figure
saveFigures <- function(plot.list, fig.folder="figures", width = 8, height = 8) {
  if (is.null(names(plot.list)))
    stop("Cannot find plot.list names !")
  for (i in 1:length(plot.list)) {
    plot <- plot.list[[i]]
    p <- gsub("p|gt", "Figure", names(plot.list)[i])
    
    if (p=="Figure2") {
      width = 12; height = 10
    } else if (p=="FigureS1" || p=="FigureS2" || p=="FigureS3") {
      width = 9; height = 9
    } else if (p=="Figure4" || p=="FigureS4") {
      width = 10; height = 12
    } else if (p=="FigureS5" || p=="FigureS6" || p=="FigureS7") {
      width = 10; height = 8
    } else if (p=="Figure5") {
      width = 7; height = 6
    } else if (p=="Figure6") {
      width = 11; height = 10
    } else if (p=="Figure7" || p=="FigureS8") {
      width = 7; height = 10
    } else if (p=="FigureS9" || p=="FigureS10" || p=="FigureS11") {
      width = 7; height = 4
    } else if (p=="FigureS12" || p=="FigureS13" || p=="FigureS14") {
      width = 7; height = 2
    } else if (p=="FigureS15") {
      width = 5; height = 2
    } else if (p=="Figure8" || p=="FigureS17" || p=="FigureS18") {
      width = 10; height = 6
    } else if (grepl("^FigureS16",p)) {
      width = 7; height = 6
    } else if (p=="Figure20" || p=="FigureS21" || p=="FigureS22") {
      width = 9; height = 9
    } else {
      stop("Cannot recognize plot varible ", names(plot.list)[i], " !")
    }
    
    if (is(plot, "gtable"))
      pdf.gtable(plot, file.path(fig.folder, paste0(p, ".pdf")), 
                 width = width, height = height)
    else if (is(plot, "ggplot"))
      pdf.ggplot(plot, file.path(fig.folder, paste0(p, ".pdf")), 
                 width = width, height = height)
    else 
      stop("Cannot find gtable/ggplot object ", class(plot), " !\nTry pdf.plot")
  }
}

replaceColNames <- function(df, pattern="\\..*", replacement="") {
  if (nchar(pattern) > 0)
    colnames(df) <- gsub(pattern, replacement, toupper(colnames(df))) 
  return(df)
}

# create pdf for plotCorrelations which cannot save to list
# Figure S19 S20 S21 S22 S23
pdfAllCorrelationsRanks <- function(ranks.list, fig.folder="figures", pattern="\\..*", 
                                    replacement="", width = 4, height = 4) {
  if (is.null(names(ranks.list)))
    stop("Cannot find ranks.list names !")
  for (i in 1:length(ranks.list)) {
    ranks <- replaceColNames(ranks.list[[i]], pattern=pattern, replacement=replacement)
    p <- gsub("p|gt", "Figure", names(ranks.list)[i])
    fig.file <- file.path(fig.folder, paste0(p, ".pdf"))
    pdf.plot(ComMA::plotCorrelations(ranks), fig.file, width = width, height = height)
    cat("Save pdf", i, " to ", fig.file, "\n")
  }
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







