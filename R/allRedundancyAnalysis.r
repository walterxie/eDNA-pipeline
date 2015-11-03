# A Caution Regarding Rules of Thumb for Variance Inflation Factors
# http://link.springer.com/article/10.1007/s11135-006-9018-6

library(vegan)
library(ggplot2)
library(grid)
library(VIF)
library(xtable)
library(scales)

if(!exists("figDir")) stop("figure folder name is missing !")
if(!exists("matrixNames")) stop("matrix names are missing !")

if(!exists("verbose")) verbose = TRUE
if(!exists("rmSingleton")) rmSingleton = TRUE
if(!exists("otuThr")) otuThr = 97
if(!exists("diss.fun")) diss.fun="beta1-1"
if(!exists("taxa.group")) taxa.group="assigned"
if(!exists("isPlot")) isPlot = FALSE # by subplot

n <- length(matrixNames) 

source("Modules/init.R", local=TRUE)
source("Modules/vif_function.R", local=TRUE)

# Prepare environmental data for analysis -------------------------------------
env <- getSampleMetaData(isPlot)
# remove columns: Plot, Elev.group
env <- env[,-c(1,9)]
env[,"ForestType"] <- gsub(":.*", "", env[,"ForestType"], ignore.case = T)
env[env=="x"] <- NA
# CM30b51, CM30b58 are removed
env <- na.omit(env)
# remove discrete values
env <- env[,-c(4,7)]
# for log transformation
env[env==0] <- 0.000000001

# Inspect data
#plot(env, gap = 0, panel = panel.smooth)
# Log transform chem variables according inspection
env[,6:13] <- log(env[,6:13])

# Inspect data again; Note that certain variables are highly correlated (e.g. EC/Organic.C/Total.N, 
# Elevation/Temperature). Also, two datasets are included - pilot vs full LBI study).
pdf(paste(workingPath, figDir, "/environmental-data.pdf", sep = ""), width=12, height=12)
plot(env, gap = 0, panel = panel.smooth, upper.panel=NULL) 
invisible(dev.off()) 

colnames(env)

rdaReducedList <- list()
rdaForwardList <- list()
rdaBackwardList <- list()
rdaDataList <- list()

minOTUs <- 200
cat("RDA Config: isPlot =", isPlot, ", rmSingleton =", rmSingleton, ", minOTUs =", minOTUs, ". \n")

count = 1
for (taxa.group in taxa.groups) {
  # Prepare OTUtable for analysis ------------------------------------------------
  for (expId in 1:(n-1)) {	
    
    communityMatrix <- getCommunityMatrixT(expId, isPlot, rmSingleton, taxa.group, minRow=minOTUs)
    
    if (is.null(communityMatrix)) {
      cat("\nSkip", taxa.group, "subset from", matrixNames[expId], ".\n") 
      next
    }
    
    # merge to match rownames
    ncol.cm <- ncol(communityMatrix)
    cm_env <- merge(communityMatrix, env, by = "row.names")
    # move 1st col Row.names to row.names
    rownames(cm_env) <- cm_env[,"Row.names"]
    cm_env <- cm_env[,-1]
    
    communityMatrix <- cm_env[,1:ncol.cm]
    cm_env <- as.data.frame(data.matrix(cm_env[,-(1:ncol.cm)]))
    
    # remove 0 row/column after merge
    communityMatrix <- prepCommunityMatrix(communityMatrix)
    
    cat("RDA of", taxa.group, "subset from", matrixNames[expId], "having", ncol(communityMatrix), "OTUs", 
        nrow(communityMatrix), "samples. \n") 
    
    # Constrained ordination ------------------------------------------------------
    rda_table <- data.frame(row.names=c("Constrained","Unconstrained"))
    anova_table <- data.frame(row.names=colnames(cm_env))
    
    # Distance-based redundancy analysis, using capscale
    # Constrained Analysis of Principal Coordinates (CAP) is an ordination method similar to Redundancy Analysis (rda).
    # DB-RDA, empty model
    rda_0 <- capscale(communityMatrix ~ 1, cm_env, distance = "jaccard")
    
    # DB-RDA, maximal model (bad idea - only use for auto model building)
    rda_1 <- capscale(communityMatrix ~ ., cm_env, distance = "jaccard")
    if (verbose) head(summary(rda_1))
    # sp = species scores, wa = site scores, bp = biplot arrows, lc = linear constraints 
    #	plot(rda_1, display = c("wa", "bp")) # Note correlation of biplot arrows
    
    # Variance inflation factor - indicates highly correlated variables
    print(vif.cca(rda_1))
    
    rda_table$Inertia <- c(round(rda_1$CCA$tot.chi, 3), round(rda_1$CA$tot.chi, 3))
    rda_table$Proportion <- c(rda_1$CCA$tot.chi/rda_1$tot.chi, rda_1$CA$tot.chi/rda_1$tot.chi)
    rda_table$Proportion <- percent(rda_table$Proportion) # %
    
    constrained_inertia <- c()
    constrained_proportion <- c()
    # Test each variable individually
    for (i in 1:length(colnames(cm_env))) {
      rda_individual <- capscale(formula = as.formula(paste("communityMatrix", colnames(cm_env)[i], sep=" ~ ")), 
                                 cm_env, distance = "jaccard")
      constrained_inertia <- c(constrained_inertia, rda_individual$CCA$tot.chi)
      constrained_proportion <- c(constrained_proportion, rda_individual$CCA$tot.chi/rda_individual$tot.chi)
    }
    
    anova_table$Inertia <- constrained_inertia
    anova_table$Proportion <- percent(constrained_proportion) # %
    
    #Compute all the single terms in the scope argument that can be added to or dropped from the model, 
    #fit those models and compute a table of the changes in fit.
    add_1 <- add1(rda_0, scope=formula(rda_1), test="perm")
    anova_table$Pr <- add_1$Pr[-1] # 1st row is <none>
    colnames(anova_table)[3] <- "Pr($>$F)"
    
    # Build model after stepwise removal of collinear variables (vif >= 10; requires vif_function.R) 
    # variance inflation factor (VIF) quantifies the severity of multicollinearity in an ordinary least squares regression analysis. 
    env_reduced <- vif_func(in_frame = cm_env)
    print(env_reduced) # Remaining variables
    
    # Build model automatically from reduced variable set
    # (Unsure how to pass env_reduced variables to capscale formula; paste() doesn't work...)
    #	rda_reduced <- capscale(communityMatrix ~ slope.degree + Mean.Temp + Northness + Eastness + 
    #							pH + C.N.ratio + NO3.N + NH4.N + Olsen.P, cm_env, distance = "jaccard")
    rda_reduced <- capscale(formula = as.formula(paste("communityMatrix", paste(env_reduced, collapse=" + "), sep=" ~ ")), 
                            cm_env, distance = "jaccard")
    if (verbose) head(summary(rda_reduced))
    anova_reduced <- anova(rda_reduced, by = "terms")
    
    rda_table$Inertia.R <- c(round(rda_reduced$CCA$tot.chi, 3), round(rda_reduced$CA$tot.chi, 3))
    rda_table$Proportion.R <- c(rda_reduced$CCA$tot.chi/rda_reduced$tot.chi, rda_reduced$CA$tot.chi/rda_reduced$tot.chi)
    rda_table$Proportion.R <- percent(rda_table$Proportion.R) # %
    
    anova_table$Reduced <- is.element(rownames(anova_table), rownames(anova_reduced))
    anova_table$Reduced[which(anova_table$Reduced==T)] <- anova_reduced$Pr[-length(anova_reduced$Pr)]
    colnames(anova_table)[4] <- "Reduced Pr($>$F)"
    
    # Choose a Model by Permutation Tests in Constrained Ordination using forward model selection
    rda_reduced_f <- ordistep(rda_0, scope = formula(rda_reduced), direction = "forward", permutations = 3999)
    
    rda_forward <- capscale(formula = as.formula(rda_reduced_f$call), data = cm_env, distance = "jaccard")
    if (verbose) head(summary(rda_forward))
    anova_forward <- anova(rda_forward, by = "terms")
    
    rda_table$Inertia.F <- c(round(rda_forward$CCA$tot.chi, 3), round(rda_forward$CA$tot.chi, 3))
    rda_table$Proportion.F <- c(rda_forward$CCA$tot.chi/rda_forward$tot.chi, rda_forward$CA$tot.chi/rda_forward$tot.chi)
    rda_table$Proportion.F <- percent(rda_table$Proportion.F) # %
    
    anova_table$Forward <- is.element(rownames(anova_table), rownames(anova_forward))
    anova_table$Forward[which(anova_table$Forward==T)] <- anova_forward$Pr[-length(anova_forward$Pr)]
    colnames(anova_table)[5] <- "Forward Pr($>$F)"
    
    # Choose a Model by Permutation Tests in Constrained Ordination using backward model selection
    rda_reduced_b <- ordistep(rda_reduced, scope = formula(rda_0), direction = "backward", permutations = 3999)
    
    rda_backward <- capscale(formula = as.formula(rda_reduced_b$call), data = cm_env, distance = "jaccard")
    if (verbose) head(summary(rda_backward))
    anova_backward <- anova(rda_backward, by = "terms")
    
    rda_table$Inertia.B <- c(round(rda_backward$CCA$tot.chi, 3), round(rda_backward$CA$tot.chi, 3))
    rda_table$Proportion.B <- c(rda_backward$CCA$tot.chi/rda_backward$tot.chi, rda_backward$CA$tot.chi/rda_backward$tot.chi)
    rda_table$Proportion.B <- percent(rda_table$Proportion.B) # %
    
    anova_table$Backward <- is.element(rownames(anova_table), rownames(anova_backward))
    anova_table$Backward[which(anova_table$Backward==T)] <- anova_backward$Pr[-length(anova_backward$Pr)]
    colnames(anova_table)[6] <- "Backward Pr($>$F)"
    
    anova_table[anova_table == 0] <- ""
    anova_table$Proportion <- gsub("%", "\\\\%", anova_table$Proportion)
    
    print(xtable(anova_table, caption = paste("Distance-based redundancy analysis and their ANOVA tests 
			in each step for the eDNA biodiversity data sets", matrixNames[expId], taxa.group), 
          label = paste("tab:rdaAnova", matrixNames[expId], taxa.group, sep = ":"), caption.placement = "top"), 
          sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)
    
    rda_table$Proportion <- gsub("%", "\\\\%", rda_table$Proportion)
    rda_table$Proportion.R <- gsub("%", "\\\\%", rda_table$Proportion.R)
    rda_table$Proportion.F <- gsub("%", "\\\\%", rda_table$Proportion.F)
    rda_table$Proportion.B <- gsub("%", "\\\\%", rda_table$Proportion.B)
    
    print(xtable(rda_table, caption = paste("The constrained and unconstrained inertia changes during 
          distance-based redundancy analysis for the eDNA biodiversity data sets", matrixNames[expId], taxa.group), 
          label = paste("tab:rda", matrixNames[expId], taxa.group, sep = ":"), caption.placement = "top"), 
          sanitize.text.function = function(x){x}, file=tableFile, append=TRUE)
    
    rdaReducedList[[ count ]] <- rda_reduced
    rdaForwardList[[ count ]] <- rda_forward
    rdaBackwardList[[ count ]] <- rda_backward
    rdaDataList[[ count ]] <- paste(matrixNames[expId], taxa.group, sep = "-")
    
    count = count + 1
    
    # figures
    fname <- paste("rda", matrixNames[expId], postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-")
    pdf(paste(workingPath, figDir, "/", fname, ".pdf", sep = ""), width=9, height=3)
    attach(mtcars)
    par(mfrow=c(1,3))
    par(mar=c(4,4,3,2), cex=0.6)
    
    rownames(rda_reduced$CCA$wa) <- substrRight(rownames(rda_reduced$CCA$wa), 4)
    rownames(rda_forward$CCA$wa) <- substrRight(rownames(rda_forward$CCA$wa), 4)
    rownames(rda_backward$CCA$wa) <- substrRight(rownames(rda_backward$CCA$wa), 4)
    
    plot(rda_reduced, display = c("wa", "bp"), main="Reduced (VIF >= 10)")
    plot(rda_forward, display = c("wa", "bp"), ylab="", main="Forward")
    plot(rda_backward, display = c("wa", "bp"), ylab="", main="Backward")
    
    invisible(dev.off())
  } # END for expId
} # END for taxa.group 

