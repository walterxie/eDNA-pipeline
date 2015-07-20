
# remove ALL objects 
rm(list=ls()) 

# change config below
sourcePath <<- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/"
setwd(sourcePath)
workingPath <<- "~/Projects/NZGO/pilot2010/pipeline/"

cat("\nConfig : set R source path", sourcePath, "\n")
cat("\nConfig : set working path", workingPath, "\n")

# print detail if TRUE
verbose <<- FALSE 

# remove all singleton from community matrix if TRUE
rmSingleton <<- FALSE 
#rmSingleton <<- TRUE  

if (rmSingleton) {
	figDir <<- "figures1"
	# write tables to a file in workingPath
	tableFile <<- paste(workingPath, "report-no-singleton.tex", sep="")
	tit <<- "remove all singletons in community matrix"	
} else {
	figDir <<- "figures"
	# write tables to a file in workingPath
	tableFile <<- paste(workingPath, "report.tex", sep="")
	tit <<- "keep all singletons in community matrix"
}
cat("\nConfig :", tit, "!\n")

# create folder for figures in workingPath
if (!file.exists(paste(workingPath, figDir, sep="")))    
    dir.create(paste(workingPath, figDir, sep=""))    

cat("\nConfig : set figures path", paste(workingPath, figDir, sep=""), "\n\n")

matrixNames <<-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name
matrixNamesNo454 <<-  c("seedlings","trees","invertebrates","birds") # need expId correct for "birds","seedlings" 

levels = rep(c("gamma","alpha","beta"),3)
qs = rep(0:2,each=3)

######## set up report latex #######
cat("\\documentclass{article}\n\n", file=tableFile, append=FALSE)
# add packages here
cat("\\usepackage[utf8]{inputenc}","\\usepackage{graphicx}","\\usepackage{caption}","\n", file=tableFile, append=TRUE, sep = "\n")
cat("\\title{Report to", tit, "}\n\n", file=tableFile, append=TRUE)
cat("\\date{\\today}","\\begin{document}", "\\maketitle", file=tableFile, append=TRUE, sep = "\n\n")

######## figures and tables #######

source("allSampleCount.r", local=TRUE)

source("allStatistics.r", local=TRUE)

if (!rmSingleton)
source("allTaxonomyPhylum.r", local=TRUE)

#source("createAllDiversitiesOTUsTable.r", local=TRUE)
#source("createAllRarefactionTable.r", local=TRUE) # only for otuThr = 97

source("allDiversitiesOTUs.r", local=TRUE)

source("allRarefactions.r", local=TRUE)

source("allWithinBetweenPlots.r", local=TRUE)

source("allElevationDiversitiesByPlots.r", local=TRUE)

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
