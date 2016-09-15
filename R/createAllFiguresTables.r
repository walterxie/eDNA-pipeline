
# remove ALL objects 
rm(list=ls()) 

# change config below
setwd("~/WorkSpace/eDNA-pipeline")

tableFile <- file.path("report.tex")

######## initialise #######
source("R/init.r", local=TRUE)

# create folder for figures 
mkdir(file.path("figures"))    

######## set up report latex #######
cat("\\documentclass{article}\n\n", file=tableFile, append=FALSE)
# add packages here
cat("\\usepackage[utf8]{inputenc}","\\usepackage{graphicx}","\\usepackage{caption}","\n", 
    file=tableFile, append=TRUE, sep = "\n")
cat("\\title{Multi-gene meta-barcoding analysis of terrestrial biodiversity on a forested island}\n\n", 
    file=tableFile, append=TRUE)
cat("\\date{\\today}","\\begin{document}", "\\maketitle", file=tableFile, append=TRUE, sep = "\n\n")

######## figures and tables #######
source("R/allStatistics.r", local=TRUE)
otu.stats <- getOTUStatistics()

source("R/allStatisticsTaxaGroup.r", local=TRUE)
tg.stats <- getTaxaGroupStatistics()

source("R/allDissimVsDistances.r", local=TRUE)
# all.dist.subplot <- getDissimVsDistances()
load("data/all.dist.subplot.RData")
# all.dist.list
# TODO UniFrac(phy, weighted = TRUE, normalized = TRUE) < 0.1
plotDistanceCorrelation(all.dist.list[["within"]])
# TODO 
plotDistanceCorrelation(all.dist.list[["elev.diff"]])
# TODO 
plotWithinBetween(all.dist.list[["within.between"]])






######## complete report latex #######

cat("\n\n\\end{document}", file=tableFile, append=TRUE, sep = "\n")

Sys.setenv(PATH="/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/texbin")
setwd(workingPath)
system(paste("pdflatex", tableFile)) 
setwd(sourcePath)

cat(paste("\n\nComplete report : ", tableFile))
