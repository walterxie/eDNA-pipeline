library(vegan)
library(Hmisc)


if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("taxa.group")) taxa.group="assigned"
if(!exists("isPlot")) isPlot = FALSE # by subplot
if(!exists("levels")) levels = rep(c("gamma","alpha","beta"),3)
if(!exists("qs")) qs = rep(0:2,each=3)

source("Modules/init.r")
source("Modules/RarefactionTable.R")

cat("Intermediate data: create rarefaction diversity table using rmSingleton =", rmSingleton, 
    ", isPlot =", isPlot, ", otuThr =", otuThr, "\n") 

n <- length(matrixNames)

filePath <- file.path(workingPath, "data", "rf")
mkdir(filePath) 

# main
for (expId in 1:(n-1)) {	
	# create rarefaction table file which is calculated by RarefactionTable.R	
  pathFileStem <- file.path(filePath, paste(matrixNames[expId], 
                            postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-"))
  
	createRDPerSampleTable(expId, isPlot, rmSingleton, taxa.group, pathFileStem, reps=10, min.size=100, verbose=verbose)
} 



