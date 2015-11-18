# A Caution Regarding Rules of Thumb for Variance Inflation Factors
# http://link.springer.com/article/10.1007/s11135-006-9018-6

library(vegan)
library(ggplot2)
library(grid)
library(VIF)
library(xtable)
library(scales)

source("Modules/init.R", local=TRUE)
source("Modules/vif_function.R", local=TRUE)

# cm.or.dist: rows are simples
# env: rows are simples, and must be same as rownames(cm.or.dist) inlcuding order
# tableFile: latex file
# verbose default TRUE
proceedRDA <- function(cm.or.dist, env, tableFile=NULL, verbose=TRUE) {
  if ( all( tolower(rownames(env)) != tolower(rownames(as.matrix(cm.or.dist))) ) ) 
    stop("Site names in community matrix and environmental data file not matched !")
  
  # Constrained ordination ------------------------------------------------------
  rda_table <- data.frame(row.names=c("Constrained","Unconstrained"))
  anova_table <- data.frame(row.names=colnames(env))
  
  # Distance-based redundancy analysis, using capscale
  # Constrained Analysis of Principal Coordinates (CAP) is an ordination method similar to Redundancy Analysis (rda).
  # DB-RDA, empty model
  rda_0 <- capscale(cm.or.dist ~ 1, env, distance = "jaccard")
  
  # DB-RDA, maximal model (bad idea - only use for auto model building)
  rda_1 <- capscale(cm.or.dist ~ ., env, distance = "jaccard")
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
  for (i in 1:length(colnames(env))) {
    rda_individual <- capscale(formula = as.formula(paste("cm.or.dist", colnames(env)[i], sep=" ~ ")), 
                               env, distance = "jaccard")
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
  env_reduced <- vif_func(in_frame = env)
  print(env_reduced) # Remaining variables
  
  # Build model automatically from reduced variable set
  # (Unsure how to pass env_reduced variables to capscale formula; paste() doesn't work...)
  #	rda_reduced <- capscale(cm.or.dist ~ slope.degree + Mean.Temp + Northness + Eastness + 
  #							pH + C.N.ratio + NO3.N + NH4.N + Olsen.P, env, distance = "jaccard")
  rda_reduced <- capscale(formula = as.formula(paste("cm.or.dist", paste(env_reduced, collapse=" + "), sep=" ~ ")), 
                          env, distance = "jaccard")
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
  
  rda_forward <- capscale(formula = as.formula(rda_reduced_f$call), data = env, distance = "jaccard")
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
  
  rda_backward <- capscale(formula = as.formula(rda_reduced_b$call), data = env, distance = "jaccard")
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
  
  printXTable(anova_table, caption = paste("Distance-based redundancy analysis and their ANOVA tests 
          in each step for the eDNA biodiversity data sets", matrixNames[expId], taxa.group), 
          label = paste("tab:rdaAnova", matrixNames[expId], taxa.group, sep = ":"), file=tableFile)
  
  rda_table$Proportion <- gsub("%", "\\\\%", rda_table$Proportion)
  rda_table$Proportion.R <- gsub("%", "\\\\%", rda_table$Proportion.R)
  rda_table$Proportion.F <- gsub("%", "\\\\%", rda_table$Proportion.F)
  rda_table$Proportion.B <- gsub("%", "\\\\%", rda_table$Proportion.B)
  
  printXTable(rda_table, caption = paste("The constrained and unconstrained inertia changes during 
          distance-based redundancy analysis for the eDNA biodiversity data sets", matrixNames[expId], taxa.group), 
          label = paste("tab:rda", matrixNames[expId], taxa.group, sep = ":"), file=tableFile)
  
  # Return a list 
  list(
    rda_reduced = rda_reduced,
    rda_forward = rda_forward,
    rda_backward = rda_backward
  )
}


