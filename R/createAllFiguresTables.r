
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
otu.stats <- getOTUStatistics(input.names, file.xtable=tableFile)

source("R/allStatisticsTaxaGroup.r", local=TRUE)
tg.stats <- getTaxaGroupStatistics(input.names, file.xtable=tableFile)
#tg.stats$phyla.list$OTUs

# use phyloseq 1.10.0, new version weighted UniFrac < 0.1
source("R/allDissimVsDistances.r", local=TRUE)
# all.dist.subplot <- getDissimVsDistances(input.names, init=FALSE, save.rdata=TRUE)

# subplot
load("data/all.dist.subplot.RData")
# all.dist.list
plotDistanceCorrelation(all.dist.list[["within"]])
plotWithinBetween(all.dist.list[["within.between"]])

# plot 
load("data/all.dist.subplot.RData")
plotDistanceCorrelation(all.dist.list[["elev.diff"]])

source("R/allTaxonomyAnalyses.r", local=TRUE)
all.counts.sums <- getAllCountsSums(init=FALSE)


######## complete report latex #######

cat("\n\n\\end{document}", file=tableFile, append=TRUE, sep = "\n")

Sys.setenv(PATH="/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/texbin")
setwd(workingPath)
system(paste("pdflatex", tableFile)) 
setwd(sourcePath)

cat(paste("\n\nComplete report : ", tableFile))
