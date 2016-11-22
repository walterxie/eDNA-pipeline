
# remove ALL objects 
rm(list=ls()) 

# change config below
setwd("~/WorkSpace/eDNA-pipeline")

######## initialise #######
source("R/init.r", local=TRUE)

# tableFile <- NULL # print to console
tableFile <- file.path("report.tex")
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
# Table 2
otu.stats <- getOTUStatistics(input.names, file.xtable=tableFile)

source("R/allStatisticsTaxaGroup.r", local=TRUE)
# Table 2
tg.stats <- getTaxaGroupStatistics(input.names, file.xtable=tableFile)
#tg.stats$phyla.list$OTUs

source("R/allTaxonomyAnalyses.r", local=TRUE)
# Figure 2
all.counts.sums <- getAllCountsSums(input.names)

# use phyloseq 1.10.0, new version weighted UniFrac < 0.1
source("R/allDissimVsDistances.r", local=TRUE)
# all.dist.subplot <- getDissimVsDistances(input.names, init=FALSE, save.rdata=TRUE)

# all.dist.list$by.plot == F 
load("data/all.dist.subplot.RData")
# Figure 3a, S1
plotDistanceCorrelation(all.dist.list[["within"]])
# Figure 3b, S2
plotDistanceCorrelation(all.dist.list[["elev.diff"]])
# Figure S3
plotWithinBetween(all.dist.list[["within.between"]])






######## complete report latex #######

cat("\n\n\\end{document}", file=tableFile, append=TRUE, sep = "\n")

Sys.setenv(PATH="/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/texbin")
setwd(workingPath)
system(paste("pdflatex", tableFile)) 
setwd(sourcePath)

cat(paste("\n\nComplete report : ", tableFile))
