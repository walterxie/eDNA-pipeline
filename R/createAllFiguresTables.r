
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
cat("\\title{Spatial patterns of soil biodiversity revealed by multi-gene meta-barcoding analysis in a forested island ecosystem}\n\n", 
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

# use phyloseq 1.10.0, new version bug weighted UniFrac < 0.1
source("R/allDissimVsDistances.r", local=TRUE)
# all.dist.subplot <- getDissimVsDistances(input.names, by.plot=FALSE, save.rdata=TRUE)
# all.dist.list$by.plot == F 
load("data/all.dist.subplot.RData")
# Figure S1, 3a is subset of S1
pS1 <- plotDistanceCorrelation(all.dist.list[["within"]])
# Figure S2, 3b is subset of S2
pS2 <- plotDistanceCorrelation(all.dist.list[["elev.diff"]])
# Figure S3
pS3 <- plotWithinBetween(all.dist.list[["within.between"]])
saveFigures(list(p2=p2, pS1=pS1, pS2=pS2, pS3=pS3))

# NMDS
source("R/allNMDS.r", local=TRUE)
# Figure 4
nmds <- getNMDS(input.names)
gt4 <- ComMA::grid_arrange_shared_legend(nmds$plot.list, input.list=T, nrow=3, 
                                         legend.position="right", widths=c(0.8, 0.2))
# Figure S4: all OTUs, Figure S5, S6, S7
genes.taxa.list <- list(gtS4=list(list("16S","all"),list("18S","all"),list("26S","all"),
                                  list("ITS","all"),list("ShCO1","all"),list("FolCO1","all")),
                        gtS5=list(list("18S","fungi"),list("26S","fungi"),
                                  list("ITS","fungi"),list("ShCO1","fungi")),
                        gtS6=list(list("18S","protists"),list("26S","protists"),
                                  list("ShCO1","protists"),list("FolCO1","protists")),
                        gtS7=list(list("18S","animals"),list("26S","animals"),
                                  list("ShCO1","animals"),list("FolCO1","animals")) )
nmds.list <- plotAllNMDS(input.names, genes.taxa.list)
nmds.list[["gt4"]] <- gt4
saveFigures(nmds.list)

# community comparison
source("R/allGeneCorrolation.r", local=TRUE)
# Mantel test (lower triangle) and Procrustes test (upper triangle) 
corrs <- getMantelAndProcrustes(input.names)
# Figure 5 heatmap: plots$heatmap; plots$mantel.mds; plots$prot.mds
plots <- plotMantelAndProcrustes(corrs)
printMantelAndProcrustes(corrs)

gene.levels=c("16S bacteria","18S animals","18S fungi","18S protists","26S animals","26S fungi","26S protists",
              "ITS fungi","COI-300 animals","COI-300 fungi","COI-300 protists","COI-650 animals","COI-650 protists")
genes.taxa=list(list("16S","bacteria"),list("18S","animals"),list("18S","fungi"),list("18S","protists"), 
                list("26S","animals"),list("26S","fungi"),list("26S","protists"),list("ITS","fungi"), 
                list("ShCO1","animals"),list("ShCO1","fungi"),list("ShCO1","protists"),
                list("FolCO1","animals"),list("FolCO1","protists"))
corrs2 <- getMantelAndProcrustes(input.names, genes.taxa=genes.taxa, order.by=gene.levels)
# Figure 6 heatmap: plots2$heatmap; plots2$mantel.mds; plots2$prot.mds
plots2 <- plotMantelAndProcrustes(corrs2, gene.levels=gene.levels)

# env data by subplots
env.subplot <- getEnvData(by.plot=F)
# Figure S8
pS8 <- ComMA::plotProcrustes(corrs$procrustes$proc.list, env.subplot, proc.list.pairs=corrs$procrustes$pairs, colour.id="Elevation")
ComMA::grid_arrange_shared_legend(pS8[[1]], input.list=T, ncol=3, nrow=5, legend.position="right", widths=c(1, 0.1, 0.1))
# Figure 7: p.all$gt7; Figure S9-15: p.all$gtS9 until gtS15
p.all <- plotAllProcrustes(corrs2$procrustes, env.subplot)

saveFigures(list(p5=plots$heatmap, p6=plots2$heatmap, pS8=pS8, gt7=p.all$gt7, gtS9=p.all$gtS9))

# ranking of sample plots
source("R/allPlotPrioritisation.r", local=TRUE)
diversities=c("gamma1","beta1","pd.alpha","sp.rich")
pp.df.list <- prioriPlotByDiversities(input.names, diversities=diversities)
#save(pp.df.list, file = "data/plot.prior.prok.euk.RData" )
# env data by plots
env.plot <- getEnvData(by.plot=T)
# Figure S16: gtS16a$heatmap
gtS16a <- ComMA::plotPrioritisation.Attribute(pp.df.list[["rank"]][["pd.alpha"]], env.plot, x.lab="Sample plot", 
                                              y.lab="Amplicon dataset", grid.widths = c(6, 3))
gtS16b <- ComMA::plotPrioritisation.Attribute(pp.df.list[["rank"]][["sp.rich"]], env.plot, x.lab="Sample plot", 
                                              y.lab="Amplicon dataset", grid.widths = c(6, 3))
gtS16c <- ComMA::plotPrioritisation.Attribute(pp.df.list[["rank"]][["gamma1"]], env.plot, x.lab="Sample plot", 
                                              y.lab="Amplicon dataset", grid.widths = c(6, 3))

pp.df2.list <- prioriPlotByDiversities(input.names, diversities=diversities, genes.taxa=genes.taxa)
#save(pp.df2.list, file = "data/plot.prior.taxa.subsets.RData" )
# Figure 8
gt8 <- ComMA::plotPrioritisation.Attribute(pp.df2.list[["rank"]][["pd.alpha"]], env.plot, x.lab="Sample plot", 
                                           y.lab="Group", grid.widths = c(8, 3))
# Figure S17
gtS17 <- ComMA::plotPrioritisation.Attribute(pp.df2.list[["rank"]][["sp.rich"]], env.plot, x.lab="Sample plot", 
                                             y.lab="Group", grid.widths = c(8, 3))
# Figure S18
gtS18 <- ComMA::plotPrioritisation.Attribute(pp.df2.list[["rank"]][["gamma1"]], env.plot, x.lab="Sample plot", 
                                             y.lab="Group", grid.widths = c(8, 3))

saveFigures(list(gtS16a=gtS16a$heatmap, gtS16b=gtS16b$heatmap, gtS16c=gtS16c$heatmap, 
                 gt8=gt8$heatmap, gtS17=gtS17$heatmap, gtS18=gtS18$heatmap))

# RDA
source("R/allRedundancyAnalysis.r", local=TRUE)
rda <- getRDA(input.names)




######## complete report latex #######

cat("\n\n\\end{document}", file=tableFile, append=TRUE, sep = "\n")

Sys.setenv(PATH="/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/texbin")
setwd(workingPath)
system(paste("pdflatex", tableFile)) 
setwd(sourcePath)

cat(paste("\n\nComplete report : ", tableFile))
