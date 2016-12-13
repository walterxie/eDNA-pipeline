
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
all.co.su <- getAllCountsSums(input.names)
printAllCountsSums(all.co.su, file.xtable=tableFile)

# use phyloseq 1.10.0, new version bug weighted UniFrac < 0.1
source("R/allDissimVsDistances.r", local=TRUE)
all.dist.list <- getDissimVsDistances(input.names, by.plot=F, save.rdata=F)
# Use subplots for Figure S1 and S3, but plots for S2 (all.dist.list$by.plot == F) 
#load("data/all.dist.subplot.RData")
# Figure S1, Figure 3a is subset of Figure S1
pS1 <- plotDistanceCorrelation(all.dist.list[["within"]])
# Figure S3
pS3 <- plotWithinBetween(all.dist.list[["within.between"]])
#load("data/all.dist.plot.RData")
all.dist.list <- getDissimVsDistances(input.names, by.plot=T, save.rdata=F)
# Figure S2, Figure 3b is subset of Figure S2
pS2 <- plotDistanceCorrelation(all.dist.list[["elev.diff"]])
saveFigures(list(p2=all.co.su$ggplot, pS1=pS1, pS2=pS2, pS3=pS3))

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
pS8 <- plotProcrustes.allOTUs(corrs$procrustes, env.subplot)
# Figure 7: p.all$gt7; Figure S9-15: p.all$gtS9 until gtS15
p.all <- plotAllProcrustes(corrs2$procrustes, env.subplot)

saveFigures(list(p5=plots$heatmap, p6=plots2$heatmap, pS8=pS8, gt7=p.all$gt7, 
                 gtS9=p.all$gtS9, gtS10=p.all$gtS10, gtS11=p.all$gtS11, gtS12=p.all$gtS12, 
                 gtS13=p.all$gtS13, gtS14=p.all$gtS14, gtS15=p.all$gtS15))

# ranking of sample plots
source("R/allPlotPrioritisation.r", local=TRUE)
diversities=c("gamma1","beta1","pd.alpha","sp.rich")
pp.df.list <- prioriPlotByDiversities(input.names, diversities=diversities)
#save(pp.df.list, file = "data/plot.prior.prok.euk.RData" )
pp.df2.list <- prioriPlotByDiversities(input.names, diversities=diversities, genes.taxa=genes.taxa)
#save(pp.df2.list, file = "data/plot.prior.taxa.subsets.RData" )
# env data by plots
env.plot <- getEnvData(by.plot=T)
# Figure S16 a b c
hm.list <- plotAllHeatmaps(list(gtS16a=pp.df.list[["rank"]][["pd.alpha"]], gtS16b=pp.df.list[["rank"]][["sp.rich"]], 
                                gtS16c=pp.df.list[["rank"]][["gamma1"]]), env.plot)
saveFigures(hm.list)
# Figure 8, S17, S18
hm.list <- plotAllHeatmaps(list(gt8=pp.df2.list[["rank"]][["pd.alpha"]], gtS17=pp.df2.list[["rank"]][["sp.rich"]], 
                                gtS18=pp.df2.list[["rank"]][["gamma1"]]), env.plot, 
                           pattern="\\.", replacement="\n", x.lab="Group", grid.widths = c(8, 1.75) )
saveFigures(hm.list)
# Figure S19 a b c
pdfAllCorrelationsRanks(list(gtS19a=pp.df.list[["rank"]][["pd.alpha"]], 
                             gtS19b=pp.df.list[["rank"]][["sp.rich"]], 
                             gtS19c=pp.df.list[["rank"]][["gamma1"]]) )
# Figure S20 S21 S22
pdfAllCorrelationsRanks(list(gtS20=pp.df2.list[["rank"]][["pd.alpha"]], 
                             gtS21=pp.df2.list[["rank"]][["sp.rich"]], 
                             gtS22=pp.df2.list[["rank"]][["gamma1"]]), 
                        pattern="\\.", replacement="\n", width = 10, height = 10 )

# RDA
source("R/allRedundancyAnalysis.r", local=TRUE)
rda <- getRDA(input.names)

pdfAllCorrelationsRanks(list(gtS23=), pattern="", width = 10, height = 10 )


######## complete report latex #######

cat("\n\n\\end{document}", file=tableFile, append=TRUE, sep = "\n")

Sys.setenv(PATH="/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/texbin")
setwd(workingPath)
system(paste("pdflatex", tableFile)) 
setwd(sourcePath)

cat(paste("\n\nComplete report : ", tableFile))
