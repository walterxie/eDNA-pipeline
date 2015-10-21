
# remove ALL objects 
rm(list=ls()) 

# change config below
sourcePath <<- "~/WorkSpace/eDNA-pipeline/R/"
setwd(sourcePath)
workingPath <<- "~/Projects/NZGO/Miseq/"

cat("\nConfig : set R source path", sourcePath, "\n")
cat("\nConfig : set working path", workingPath, "\n")

# remove all singleton from community matrix if TRUE
#rmSingleton <<- FALSE 
rmSingleton <<- TRUE  

if (rmSingleton) {
	figDir <<- "figures"
	# write tables to a file in workingPath
	tableFile <<- paste(workingPath, "report.tex", sep="")
	tit <<- "remove all singletons"	
} else {
	figDir <<- "figures1"
	# write tables to a file in workingPath
	tableFile <<- paste(workingPath, "report-singleton.tex", sep="")
	tit <<- "keep all singletons"
}
cat("\nConfig :", tit, "!\n")

######## initialise #######
source("Modules/init.r", local=TRUE)

# create folder for figures in workingPath
mkdir(file.path(workingPath, figDir))    

######## set up analysis #######
matrixNames <<-  c("16S", "18S", "26S", "ITS", "ShCO1", "FolCO1", "Vegetation") # for file name and dataset names in figures and tables
taxaFiles <<- c("16S_taxonomy_table.txt", "18S_taxonomy_table.txt", "26S_taxonomy_table.txt", "ITS_taxonomy_table.txt", 
				"ShCO1_taxonomy_table.txt", "FolCO1_taxonomy_table.txt") # only eDNA

otuThr = 97
levels = rep(c("gamma","alpha","beta"),3)
qs = rep(0:2,each=3)

verbose <<- TRUE # only print "Upload community matrix" 1st time 

######## set up report latex #######
cat("\\documentclass{article}\n\n", file=tableFile, append=FALSE)
# add packages here
cat("\\usepackage[utf8]{inputenc}","\\usepackage{graphicx}","\\usepackage{caption}","\n", file=tableFile, append=TRUE, sep = "\n")
cat("\\title{eDNA data pipeline report to", tit, "}\n\n", file=tableFile, append=TRUE)
cat("\\date{\\today}","\\begin{document}", "\\maketitle", file=tableFile, append=TRUE, sep = "\n\n")

######## figures and tables #######
# print detail if TRUE
isPlot <<- FALSE # by subplot
source("allSampleCount.r", local=TRUE)

verbose <<- FALSE
isPlot <<- TRUE # by plot
source("allStatistics.r", local=TRUE)

source("allTaxonomyPhylum.r", local=TRUE)

isPlot <<- FALSE
# create all rarefaction table
#source("createAllRarefactionTable.r", local=TRUE)
source("allRarefactions.r", local=TRUE)

# create all beta1-1 matrix for the rest of analyses 
#source("createAllBeta1Minus1.r", local=TRUE)

source("allWithinBetweenPlots.r", local=TRUE)

isPlot <<- TRUE
source("allElevationDiversitiesByPlots.r", local=TRUE)

isPlot <<- FALSE
source("allMDSBySubplots.r", local=TRUE)

source("allProcrustes.r", local=TRUE)

source("allGeneCorrolation.r", local=TRUE)

######## supplementary #######

source("allElevationAlpha.r", local=TRUE)

source("allMaxDivCombOfPlots.r", local=TRUE)

cat("\n\n\\clearpage\n\n", file=tableFile, append=TRUE) # too much floating table/figures

source("allMaxRemainedDiversity.r", local=TRUE)

cat("\n\n\\clearpage\n\n", file=tableFile, append=TRUE) # too much floating table/figures

source("allRedundancyAnalysis.r", local=TRUE)

source("allRedundancyAnalysisPlants.r", local=TRUE)

######## complete report latex #######

cat("\n\n\\end{document}", file=tableFile, append=TRUE, sep = "\n")

Sys.setenv(PATH="/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/texbin")
setwd(workingPath)
system(paste("pdflatex", tableFile)) 
setwd(sourcePath)

cat(paste("\n\nComplete report : ", tableFile))
