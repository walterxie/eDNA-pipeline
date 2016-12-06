
# remove ALL objects 
rm(list=ls()) 

# change config below
setwd("~/WorkSpace/eDNA-pipeline")

######## initialise #######
source("R/init.r", local=TRUE)

# tableFile <- NULL # print to console
tableFile <- file.path("report.tex")
fig.folder <- file.path("figures")
# create folder for figures 
mkdir(fig.folder)    

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
p2 <- getAllCountsSums(input.names)
ggsave(p2, file = file.path(fig.folder, "Overall_taxonomy_OTUs_reads_by_phylum.pdf"), 
       width = 340, height = 260, units = "mm")


# use phyloseq 1.10.0, new version bug weighted UniFrac < 0.1
source("R/allDissimVsDistances.r", local=TRUE)
# all.dist.subplot <- getDissimVsDistances(input.names, by.plot=FALSE, save.rdata=TRUE)

# all.dist.list$by.plot == F 
load("data/all.dist.subplot.RData")
# Figure 3a, S1
ps1 <- plotDistanceCorrelation(all.dist.list[["within"]])
# Figure 3b, S2
ps2 <- plotDistanceCorrelation(all.dist.list[["elev.diff"]])
# Figure S3
ps3 <- plotWithinBetween(all.dist.list[["within.between"]])

# NMDS
source("R/allNMDS.r", local=TRUE)
# Figure 4
nmds <- getNMDS(input.names)
ComMA::grid_arrange_shared_legend(nmds$plot.list, input.list=T, nrow=3, legend.position="right")
# Figure S4: all OTUs
nmds <- getNMDS(input.names, genes.taxa=list(list("16S","all"),list("18S","all"),list("26S","all"),
                                             list("ITS","all"),list("ShCO1","all"),list("FolCO1","all")) )
ComMA::grid_arrange_shared_legend(nmds$plot.list, input.list=T, nrow=3, legend.position="right")
# Figure S5
nmds <- getNMDS(input.names, genes.taxa=list(list("18S","fungi"),list("26S","fungi"),
                                             list("ITS","fungi"),list("ShCO1","fungi")) )
ComMA::grid_arrange_shared_legend(nmds$plot.list, input.list=T, legend.position="right")
# Figure S6
nmds <- getNMDS(input.names, genes.taxa=list(list("18S","protists"),list("26S","protists"),
                                             list("ShCO1","protists"),list("FolCO1","protists")) )
ComMA::grid_arrange_shared_legend(nmds$plot.list, input.list=T, legend.position="right")
# Figure S7
nmds <- getNMDS(input.names, genes.taxa=list(list("18S","animals"),list("26S","animals"),
                                             list("ShCO1","animals"),list("FolCO1","animals")) )
ComMA::grid_arrange_shared_legend(nmds$plot.list, input.list=T, legend.position="right")

# community comparison
source("R/allGeneCorrolation.r", local=TRUE)
# Mantel test (lower triangle) and Procrustes test (upper triangle) 
corrs <- getMantelAndProcrustes(input.names)
# Figure 5 heatmap
plots <- plotMantelAndProcrustes(corrs)
# plots$heatmap; plots$mantel.mds; plots$prot.mds
printMantelAndProcrustes(corrs)

gene.levels=c("16S bacteria","18S animals","18S fungi","18S protists","26S animals","26S fungi","26S protists",
              "ITS fungi","COI-300 animals","COI-300 fungi","COI-300 protists","COI-650 animals","COI-650 protists")
corrs2 <- getMantelAndProcrustes(input.names, 
                                genes.taxa=list(list("16S","bacteria"),list("18S","animals"),list("18S","fungi"),
                                                list("18S","protists"), list("26S","animals"),list("26S","fungi"),
                                                list("26S","protists"),list("ITS","fungi"), list("ShCO1","animals"),
                                                list("ShCO1","fungi"),list("ShCO1","protists"),
                                                list("FolCO1","animals"),list("FolCO1","protists")),
                                order.by=gene.levels)
# Figure 6 heatmap
plots <- plotMantelAndProcrustes(corrs2, gene.levels=gene.levels)

# env data by subplots
env.subplot <- getEnvData(by.plot=F)
# Figure S8
pS8 <- ComMA::plotProcrustes(corrs$procrustes$proc.list, env.subplot, proc.list.pairs=corrs$procrustes$pairs, colour.id="Elevation")
ComMA::grid_arrange_shared_legend(pS8[[1]], input.list=T, ncol=3, nrow=5, legend.position="right", widths=c(1, 0.1, 0.1))
# Figure 7, Figure S9-15
p.all <- plotAllProcrustes(corrs2$procrustes, env.subplot)
p.all$gt7 #p.all$gtS9 until gtS15





######## complete report latex #######

cat("\n\n\\end{document}", file=tableFile, append=TRUE, sep = "\n")

Sys.setenv(PATH="/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/texbin")
setwd(workingPath)
system(paste("pdflatex", tableFile)) 
setwd(sourcePath)

cat(paste("\n\nComplete report : ", tableFile))
