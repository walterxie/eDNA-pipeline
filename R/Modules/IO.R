# Author: Walter Xie
# Accessed on 21 Aug 2015

library(tools)
# create folder in path if not exist, such as setup.dir(path, "figures")
mkdir <- function(subDir.path) {
  if (!file.exists(subDir.path)) {    
    dir.create(subDir.path)    
  }
  cat("\nConfig : setup", subDir.path, "\n")
}


write.cm <- function(df, file){
  if (tolower(file_ext(file))=="csv") {
    write.csv(df, file, quote=FALSE)
  } else { # .tsv .txt
    #write.table bug: mistake to start col names from the 1st cell
    cat("",colnames(df),file=file,sep="\t")
    cat("\n",file=file, append=TRUE)
    write.table(df, file, sep ="\t", quote=FALSE, col.names=FALSE, append=TRUE)
  }
}

# hasGroup, specify if values in the last column are groups, it affects how to process matrix
# hasGroup=TRUE, return a data frame by removing last column (groups), 
# and another 1-column data frame for the last column (groups). 
# return a list, where data frame communityMatrix, cols are samples, rows are OTUs/taxa, no preprocessing
# data frame groups may be NULL depending on hasGroup
readCommunityMatrixFile <- function(file=file, hasGroup) { 
  if(missing(hasGroup)) hasGroup=FALSE

  communityMatrix <- readFile(file)
  cat("\nUpload community matrix file : ", ncol(communityMatrix), "columns,", nrow(communityMatrix), "rows, from", file, "\n") 
  groups <- NULL
  
  # set NA (empty cell) to 0
  communityMatrix[is.na(communityMatrix)] <- 0
  if(hasGroup) {
    groups <- data.frame(row.names=rownames(communityMatrix), groups=communityMatrix[,ncol(communityMatrix)])
    communityMatrix <- communityMatrix[,-ncol(communityMatrix)]
    groups.unique <- unique(groups[,1])
    cat("Split last column to data frame groups that contains", length(groups.unique), 
        "groups : ", paste(groups.unique, collapse=", "), ".\n") 
  } 

  # Return a list 
  list(
    communityMatrix = communityMatrix,
    groups = groups
  )
}

readTaxaFile <- function(file=file) { 
  taxa <- readFile(file)
  cat("\nUpload taxa file : ", ncol(taxa), "columns,", nrow(taxa), "rows, from", file, "\n") 
  return(taxa)
}

readEnvDataFile <- function(file=file) { 
  envData <- readFile(file)
  cat("\nUpload environmental data file : ", ncol(envData), "columns,", nrow(envData), "rows, from", file, "\n") 
  return(envData)
}

readFile <- function(file=file, sep) { 
  if(missing(sep)) sep="\t" # only work for non csv file
  
  if (tolower(file_ext(file))=="csv") {
    df <- read.csv(file, head=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
  } else {
    # be careful read.table bug   
    df <- read.table(file, sep=sep, header=T, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)  
  }	
  return(df)
}

printXTable <- function(df, caption, label, file=NULL, align = NULL, digits = NULL) {
  if (is.null(file)) {
    print(xtable(df, caption = caption, label = label, caption.placement = "top", 
          align = align, digits = digits),
          sanitize.text.function = function(x){x})
  } else {
    print(xtable(df, caption = caption, label = label, caption.placement = "top", 
          align = align, digits = digits),
          sanitize.text.function = function(x){x}, file=file, append=TRUE)
  }
}


