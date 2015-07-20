#### input otu table ####
# rows are otu name
# columns are subplot name

## for 454 data, if matrixName contains "by-plot" then merge nearby 2 columns into 1
## but its and 16s are special:
## its 16s has 2 data set and 20 subplot each, 40 data columns and 1 extra column for otu name
## cm_subplots has 20 columns only, otu name in rownames 

# validate input
if(!exists("inputCM")) stop("inputCM is missing !")

if (grepl("\\.csv", inputCM, ignore.case = TRUE)) {
    cm_subplots <- read.csv(file=inputCM, head=TRUE, sep=",", row.names=1, check.names=FALSE)
} else {   
	cm_subplots <- read.table(inputCM, header=T, row.names=1, check.names=FALSE)  
}

print(paste("original data: columns = ", ncol(cm_subplots), ", rows = ", nrow(cm_subplots), ", file = ", inputCM, sep=""))  

# merge replicate data set
if (replicate > 1) {
    cm <- cm_subplots
    n <- length(cm) / 2
	cm_subplots <- data.frame(row.names=rownames(cm))
	for (i in 1:n) {
		newName <- substring(colnames(cm)[i],4); # hard code for names
		cm_subplots[,newName] <- rowSums(cm[,c(i,i+n)]); # hard code for specific format
	} 
}

# 2 subplots within the plot should be next each other 
for (i in seq(1,ncol(cm_subplots),by=2)) {
    plotName1 <- unlist(strsplit(colnames(cm_subplots)[i], "-", fixed = TRUE))[1] # hard code for names Plot1-B
    plotName2 <- unlist(strsplit(colnames(cm_subplots)[i+1], "-", fixed = TRUE))[1] 
	if(plotName1 != plotName2) stop("Two subplots within the same plot should be next each other !") 
}

# if matrixName contains "by-plot" then merge nearby 2 columns into 1
if (grepl("by-plot", matrixName, ignore.case = TRUE)) {
	n <- length(cm_subplots) / 2

	communityMatrix <- data.frame(row.names=rownames(cm_subplots))
	for (i in 1:n) {
		newName <- unlist(strsplit(colnames(cm_subplots)[i*2], "-", fixed = TRUE))[1] # hard code for names 1-B
		communityMatrix[,newName] <- rowSums(cm_subplots[,c(i*2-1,i*2)]); # hard code for specific format
	}  
	
	print(paste("community matrix name contains \'by-plot\', merge matrix nearby 2 columns into 1", sep=""))
} else {
    communityMatrix <- cm_subplots    
}

if (verbose) 
print(paste("community matrix ", matrixName, ": columns = ", ncol(communityMatrix), ", rows = ", nrow(communityMatrix), sep=""))  

communityMatrix <- as.matrix(communityMatrix) 
