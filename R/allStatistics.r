

if(!exists("tableFile")) stop("table file is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("levels")) levels = rep(c("gamma","alpha","beta"),3)
if(!exists("qs")) qs = rep(0:2,each=3)

otusRowNames <- c("Reads", "OTUs", "Singleton", "Coupleton") # cannot get unique reads from CM
divRowNames <- c("$^0D_\\gamma$","$^0D_\\alpha$","$^0D_\\beta$","$^1D_\\gamma$","$^1D_\\alpha$","$^1D_\\beta$",
				"$^2D_\\gamma$","$^2D_\\alpha$","$^2D_\\beta$")
lotus = length(otusRowNames)

n <- length(matrixNames) 

#source("Modules/init.R", local=TRUE)

######## summary by datasets including singletons #######

cat("\nTable: summary by datasets including singletons: otuThr =", otuThr, "\n") 

# last column total
by.datasets <- data.frame(row.names=c(otusRowNames, divRowNames), stringsAsFactors=FALSE)

for (expId in 1:n) {	
  communityMatrix <- getCommunityMatrixT(expId, TRUE, FALSE) # always including singletons
  
  diversity_table <- diversity.df(communityMatrix)
  
  by.datasets[1,expId] <- sum(communityMatrix)
  by.datasets[2,expId] <- ncol(communityMatrix)
  by.datasets[3,expId] <- sum(colSums(communityMatrix)==1)
  by.datasets[4,expId] <- sum(colSums(communityMatrix)==2)
  for (i in 1:length(unlist(diversity_table))) {
    by.datasets[i+lotus,expId] <- unlist(diversity_table)[i]
  }
  colnames(by.datasets)[expId] <- matrixNames[expId]
}

by.datasets[1:4,"Total"] <- rowSums(by.datasets)[1:4]
#by.datasets[1:4,] <- round(by.datasets[1:4,], 0)
#by.datasets[5:nrow(by.datasets),] <- round(by.datasets[5:nrow(by.datasets),], 2)
by.datasets <- round(by.datasets, 2)
by.datasets <- format(by.datasets, big.mark=",", scientific=F)
for (i in 1:lotus)
  by.datasets[i,] <- gsub(".00", "", by.datasets[i,])

align.v <- rep("r", ncol(by.datasets) + 1)
printXTable(by.datasets, caption = paste("Table of sequence statistics for eDNA data sets", 
                                         paste(matrixNames, collapse = ", ")),
            align = align.v, label = "tab:biodiveDNA", file=tableFile)

