library(vegan)
library(ggplot2)
library(grid)
library(VIF)
library(xtable)
library(scales)

# change config below
sourcePath <- "~/svn/compevol/research/NZGenomicObservatory/Metabarcoding/R/Modules/"
setwd(sourcePath)

workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"
experiments <-  c("16S", "18S", "trnL", "ITS", "COI", "COI-spun") # only for cm file name and folder name   
matrixNames <-  experiments # 

stringBySubOrPlot <- "-by-subplot"
otuThr <- 97

source("init.R", local=TRUE)
source("vif_function.R", local=TRUE)

# print detail of analysis if TRUE
verbose <- FALSE

# Prepare environmental data for analysis -------------------------------------
envdata <- read.table(paste(workingPath, "data/LJ12027.txt", sep=""), sep="\t", header=T, row.names=1)
elev <- read.table(paste(workingPath, "data/plot_elevations.txt", sep=""), sep="\t", header=T, row.names=1)

rownames(envdata) <- gsub("CM30c30", "Plot9", rownames(envdata))
rownames(envdata) <- gsub("LB1", "Plot10", rownames(envdata))
rownames(elev) <- gsub("CM30C30", "Plot9", rownames(elev))
rownames(elev) <- gsub("LB1", "Plot10", rownames(elev))

envdata <- envdata[order(rownames(envdata)), ]
elev <- elev[order(rownames(elev)), ]

if ( all( tolower(rownames(envdata)) != tolower(rownames(elev)) ) ) 
	stop("Site names in environmental data and elevation data file not matched !")
	
elev_env <- cbind(elev, envdata[,-1])	
# Inspect data
#plot(elev_env, gap = 0, panel = panel.smooth)
# Log transform chem variables according inspection
elev_env[,c(7,9,10,11,13,14,15)] <- log(elev_env[,c(7,9,10,11,13,14,15)])

# Inspect data again; Note that certain variables are highly correlated (e.g. EC/Organic.C/Total.N, 
# Elevation/Temperature). Also, two datasets are included - pilot vs full LBI study).
pdf(paste(workingPath, "figures/environmental-data.pdf", sep = ""), width=12, height=12)
plot(elev_env, gap = 0, panel = panel.smooth) 
invisible(dev.off()) 

colnames(elev_env)[1] <- "Elevation"
colnames(elev_env)[2] <- "Slope"
colnames(elev_env)[3] <- "Aspect"
#colnames(elev_env)[4] <- "Temperature"
colnames(elev_env)

rownames(elev_env) <- gsub("Plot", "", rownames(elev_env))

rdaReducedList <- list()
rdaForwardList <- list()
rdaBackwardList <- list()

# Prepare OTUtable for analysis ------------------------------------------------
for (expId in 1:length(experiments)) {	
	communityMatrix <- init(expId, otuThr, stringBySubOrPlot)
	rownames(communityMatrix) <- gsub("CM30C30", "Plot9", rownames(communityMatrix))
	rownames(communityMatrix) <- gsub("LB1", "Plot10", rownames(communityMatrix))
	
	communityMatrix <- communityMatrix[order(rownames(communityMatrix)), ]	
	rownames(communityMatrix) <- gsub("Plot", "", rownames(communityMatrix))
	
	if ( all( tolower(rownames(elev_env)) != tolower(rownames(communityMatrix)) ) ) 
		stop("Site names in community matrix and environmental-elevation data file not matched !")

	# Constrained ordination ------------------------------------------------------
	rda_table <- data.frame(row.names=c("Constrained","Unconstrained"))
    anova_table <- data.frame(row.names=colnames(elev_env))

	# Distance-based redundancy analysis, using capscale
	# Constrained Analysis of Principal Coordinates (CAP) is an ordination method similar to Redundancy Analysis (rda).
	# DB-RDA, empty model
	rda_0 <- capscale(communityMatrix ~ 1, elev_env, distance = "jaccard")

	# DB-RDA, maximal model (bad idea - only use for auto model building)
	rda_1 <- capscale(communityMatrix ~ ., elev_env, distance = "jaccard")
	if (verbose) head(summary(rda_1))
	# sp = species scores, wa = site scores, bp = biplot arrows, lc = linear constraints 
#	plot(rda_1, display = c("wa", "bp")) # Note correlation of biplot arrows

	# Variance inflation factor - indicates highly correlated variables
	vif.cca(rda_1)

	rda_table$Inertia <- c(round(rda_1$CCA$tot.chi, 3), round(rda_1$CA$tot.chi, 3))
	rda_table$Proportion <- c(rda_1$CCA$tot.chi/rda_1$tot.chi, rda_1$CA$tot.chi/rda_1$tot.chi)
	rda_table$Proportion <- percent(rda_table$Proportion) # %
	
	constrained_inertia <- c()
	constrained_proportion <- c()
	# Test each variable individually
	for (i in 1:length(colnames(elev_env))) {
		rda_individual <- capscale(formula = as.formula(paste("communityMatrix", colnames(elev_env)[i], sep=" ~ ")), 
									elev_env, distance = "jaccard")
		constrained_inertia <- c(constrained_inertia, rda_individual$CCA$tot.chi)
		constrained_proportion <- c(constrained_proportion, rda_individual$CCA$tot.chi/rda_individual$tot.chi)
	}
			
	anova_table$Inertia <- constrained_inertia
	anova_table$Proportion <- percent(constrained_proportion) # %

	#Compute all the single terms in the scope argument that can be added to or dropped from the model, 
	#fit those models and compute a table of the changes in fit.
	add_1 <- add1(rda_0, scope=formula(rda_1), test="perm")
	anova_table$Pr <- add_1$Pr[-1] # 1st row is <none>
    colnames(anova_table)[3] <- "Pr(>F)"
    
	# Build model after stepwise removal of collinear variables (vif > 10; requires vif_function.R) 
	# variance inflation factor (VIF) quantifies the severity of multicollinearity in an ordinary least squares regression analysis. 
	env_reduced <- vif_func(in_frame = elev_env)
	env_reduced # Remaining variables

	# Build model automatically from reduced variable set
	# (Unsure how to pass env_reduced variables to capscale formula; paste() doesn't work...)
#	rda_reduced <- capscale(communityMatrix ~ slope.degree + Mean.Temp + Northness + Eastness + 
#							pH + C.N.ratio + NO3.N + NH4.N + Olsen.P, elev_env, distance = "jaccard")
	rda_reduced <- capscale(formula = as.formula(paste("communityMatrix", paste(env_reduced, collapse=" + "), sep=" ~ ")), 
							elev_env, distance = "jaccard")
	if (verbose) head(summary(rda_reduced))
	anova_reduced <- anova(rda_reduced, by = "terms")

	rda_table$Inertia.R <- c(round(rda_reduced$CCA$tot.chi, 3), round(rda_reduced$CA$tot.chi, 3))
	rda_table$Proportion.R <- c(rda_reduced$CCA$tot.chi/rda_reduced$tot.chi, rda_reduced$CA$tot.chi/rda_reduced$tot.chi)
	rda_table$Proportion.R <- percent(rda_table$Proportion.R) # %
	
	anova_table$Reduced <- is.element(rownames(anova_table), rownames(anova_reduced))
	anova_table$Reduced[which(anova_table$Reduced==T)] <- anova_reduced$Pr[-length(anova_reduced$Pr)]
	colnames(anova_table)[4] <- "Reduced Pr(>F)"

    # Choose a Model by Permutation Tests in Constrained Ordination using forward model selection
	rda_reduced_f <- ordistep(rda_0, scope = formula(rda_reduced), direction = "forward")

	rda_forward <- capscale(formula = as.formula(rda_reduced_f$call), data = elev_env, distance = "jaccard")
	if (verbose) head(summary(rda_forward))
	anova_forward <- anova(rda_forward, by = "terms")

	rda_table$Inertia.F <- c(round(rda_forward$CCA$tot.chi, 3), round(rda_forward$CA$tot.chi, 3))
	rda_table$Proportion.F <- c(rda_forward$CCA$tot.chi/rda_forward$tot.chi, rda_forward$CA$tot.chi/rda_forward$tot.chi)
	rda_table$Proportion.F <- percent(rda_table$Proportion.F) # %
	
	anova_table$Forward <- is.element(rownames(anova_table), rownames(anova_forward))
	anova_table$Forward[which(anova_table$Forward==T)] <- anova_forward$Pr[-length(anova_forward$Pr)]
	colnames(anova_table)[5] <- "Forward Pr(>F)"

	# Choose a Model by Permutation Tests in Constrained Ordination using backward model selection
	rda_reduced_b <- ordistep(rda_reduced, scope = formula(rda_0), direction = "backward")

	rda_backward <- capscale(formula = as.formula(rda_reduced_b$call), data = elev_env, distance = "jaccard")
	if (verbose) head(summary(rda_backward))
	anova_backward <- anova(rda_backward, by = "terms")

	rda_table$Inertia.B <- c(round(rda_backward$CCA$tot.chi, 3), round(rda_backward$CA$tot.chi, 3))
	rda_table$Proportion.B <- c(rda_backward$CCA$tot.chi/rda_backward$tot.chi, rda_backward$CA$tot.chi/rda_backward$tot.chi)
	rda_table$Proportion.B <- percent(rda_table$Proportion.B) # %
	
	anova_table$Backward <- is.element(rownames(anova_table), rownames(anova_backward))
	anova_table$Backward[which(anova_table$Backward==T)] <- anova_backward$Pr[-length(anova_backward$Pr)]
	colnames(anova_table)[6] <- "Backward Pr(>F)"
	
	anova_table[anova_table == 0] <- ""

	print(xtable(anova_table), sanitize.text.function = function(x){x})
	
	print(xtable(rda_table), sanitize.text.function = function(x){x})

	rdaReducedList[[ expId ]] <- rda_reduced
	rdaForwardList[[ expId ]] <- rda_forward
	rdaBackwardList[[ expId ]] <- rda_backward

	# figures
#	pdf(paste(workingPath, "figures/rda-reduced-", experiments[expId], "-", otuThr, ".pdf", sep = ""), width=6, height=6)
#	plot(rda_reduced, display = c("wa", "bp"))
#	invisible(dev.off()) 
	
#	pdf(paste(workingPath, "figures/rda-forward-", experiments[expId], "-", otuThr, ".pdf", sep = ""), width=6, height=6)
#   plot(rda_forward, display = c("wa", "bp"))
#    invisible(dev.off())

#	pdf(paste(workingPath, "figures/rda-backward-", experiments[expId], "-", otuThr, ".pdf", sep = ""), width=6, height=6)
#    plot(rda_backward, display = c("wa", "bp"))
#    invisible(dev.off())
}

subTitles <- c("(a)","(b)","(c)","(d)","(e)","(f)")

########## Reduced ###########

pdf(paste(workingPath, "figures/all-rda-reduced-", otuThr, ".pdf", sep = ""), width=6, height=9)	
attach(mtcars)
par(mfrow=c(3,2))
	
for (expId in 1:length(experiments)) {		   
	if (expId > 4) {		
		xlab="CAP1"
		par(mar=c(4,4,3,2)) 
	} else {
		xlab=""
		par(mar=c(4,4,3,2))
	}
	
	if (expId %% 2 != 0) {		
		ylab="CAP2"		
	} else {
		ylab=""
	}
	
	plot(rdaReducedList[[ expId ]], display = c("wa", "bp"), xlab=xlab, ylab=ylab, main=paste(subTitles[expId], matrixNames[expId], sep = " "))		  
}
invisible(dev.off()) 

########## Forward ###########

pdf(paste(workingPath, "figures/all-rda-forward-", otuThr, ".pdf", sep = ""), width=6, height=9)	
attach(mtcars)
par(mfrow=c(3,2))
	
for (expId in 1:length(experiments)) {		   
	if (expId > 4) {		
		xlab="CAP1"
		par(mar=c(4,4,3,2)) 
	} else {
		xlab=""
		par(mar=c(4,4,3,2))
	}
	
	if (expId %% 2 != 0) {		
		ylab="CAP2"		
	} else {
		ylab=""
	}
	
	plot(rdaForwardList[[ expId ]], display = c("wa", "bp"), xlab=xlab, ylab=ylab, main=paste(subTitles[expId], matrixNames[expId], sep = " "))		  
}
invisible(dev.off()) 

########## Backward ###########

pdf(paste(workingPath, "figures/all-rda-backward-", otuThr, ".pdf", sep = ""), width=6, height=9)	
attach(mtcars)
par(mfrow=c(3,2))
	
for (expId in 1:length(experiments)) {		   
	if (expId > 4) {		
		xlab="CAP1"
		par(mar=c(4,4,3,2)) 
	} else {
		xlab=""
		par(mar=c(4,4,3,2))
	}
	
	if (expId %% 2 != 0) {		
		ylab="CAP2"		
	} else {
		ylab=""
	}
	
	plot(rdaBackwardList[[ expId ]], display = c("wa", "bp"), xlab=xlab, ylab=ylab, main=paste(subTitles[expId], matrixNames[expId], sep = " "))		  
}
invisible(dev.off()) 
