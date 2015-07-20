#### input otu table ####
# rows are otu name
# columns are subplot name

## data column is plot, such as Sanger sequence or observation data from NZGO database

# validate input
if(!exists("inputCM")) stop("inputCM is missing !")

if (grepl("\\.csv", inputCM, ignore.case = TRUE)) {
    communityMatrix <- read.csv(file=inputCM, head=TRUE, sep=",", row.names=1, check.names=FALSE)
} else {   
	communityMatrix <- read.table(inputCM, header=T, row.names=1, check.names=FALSE)  
}

print(paste("community matrix ", matrixName, ": columns = ", ncol(communityMatrix), ", rows = ", nrow(communityMatrix), ", file = ", inputCM, sep=""))  

communityMatrix <- as.matrix(communityMatrix) 
