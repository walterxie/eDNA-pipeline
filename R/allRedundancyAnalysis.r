library(vegan)
library(ggplot2)
library(grid)
library(VIF)

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
plot(elev_env, gap = 0, panel = panel.smooth)

# Prepare OTUtable for analysis ------------------------------------------------
for (expId in 1:length(experiments)) {	
	communityMatrix <- init(expId, otuThr, stringBySubOrPlot)
	rownames(communityMatrix) <- gsub("CM30C30", "Plot9", rownames(communityMatrix))
	rownames(communityMatrix) <- gsub("LB1", "Plot10", rownames(communityMatrix))
	communityMatrix <- communityMatrix[order(rownames(communityMatrix)), ]
	
	if ( all( tolower(rownames(elev_env)) != tolower(rownames(communityMatrix)) ) ) 
		stop("Site names in community matrix and environmental-elevation data file not matched !")

	# Constrained ordination ------------------------------------------------------
	rownames(elev_env) <- gsub("Plot", "", rownames(elev_env))
	rownames(communityMatrix) <- gsub("Plot", "", rownames(communityMatrix))

	# Distance-based redundancy analysis, using capscale

	# DB-RDA, empty model
	rda_0 <- capscale(communityMatrix ~ 1, elev_env, distance = "jaccard")

	# DB-RDA, maximal model (bad idea - only use for auto model building)
	rda_1 <- capscale(communityMatrix ~ ., elev_env, distance = "jaccard")
	head(summary(rda_1))
	# sp = species scores, wa = site scores, bp = biplot arrows, lc = linear constraints 
	plot(rda_1, display = c("wa", "bp")) # Note correlation of biplot arrows

	# Variance inflation factor - indicates highly correlated variables
	vif.cca(rda_1)

	# Test each variable individually
	capscale(communityMatrix ~ Elevation.m, elev_env, distance = "jaccard")
	capscale(communityMatrix ~ slope.degree, elev_env, distance = "jaccard")
	capscale(communityMatrix ~ aspect.degree, elev_env, distance = "jaccard")
	capscale(communityMatrix ~ Mean.Temp, elev_env, distance = "jaccard")
	capscale(communityMatrix ~ Northness, elev_env, distance = "jaccard")
	capscale(communityMatrix ~ Eastness, elev_env, distance = "jaccard")
	capscale(communityMatrix ~ Water.Content, elev_env, distance = "jaccard")
	capscale(communityMatrix ~ pH, elev_env, distance = "jaccard")
	capscale(communityMatrix ~ EC, elev_env, distance = "jaccard")
	capscale(communityMatrix ~ Organic.C, elev_env, distance = "jaccard")
	capscale(communityMatrix ~ Total.N, elev_env, distance = "jaccard")
	capscale(communityMatrix ~ C.N.ratio, elev_env, distance = "jaccard")
	capscale(communityMatrix ~ NO3.N, elev_env, distance = "jaccard")
	capscale(communityMatrix ~ NH4.N, elev_env, distance = "jaccard")
	capscale(communityMatrix ~ Olsen.P, elev_env, distance = "jaccard")

	add1(rda_0, scope=formula(rda_1), test="perm")

	# Build model automatically (do forward and backward methods give similar results?)
	rda_f <- ordistep(rda_0, scope = formula(rda_1), direction = "forward")
	rda_b <- ordistep(rda_1, scope = formula(rda_0), direction = "backward")
	head(summary(rda_f))
	head(summary(rda_b))
	anova(rda_f, by = "terms")
	anova(rda_b, by = "terms")
	plot(rda_f, display = c("wa", "bp"))
	plot(rda_b, display = c("wa", "bp"))

	# Build model after stepwise removal of collinear variables (vif > 10; requires vif_function.R) 
	env_reduced <- vif_func(in_frame = elev_env)
	env_reduced # Remaining variables

	# Build model automatically from reduced variable set
	# (Unsure how to pass env_reduced variables to capscale formula; paste() doesn't work...)
	rda_reduced <- capscale(communityMatrix ~ slope.degree + aspect.degree + Mean.Temp + Northness + 
							pH + EC + C.N.ratio + NO3.N + Olsen.P, elev_env, distance = "jaccard")
	head(summary(rda_reduced))
	anova(rda_reduced, by = "terms")
	plot(rda_reduced, display = c("wa", "bp"))

	rda_reduced_f <- ordistep(rda_0, scope = formula(rda_reduced), direction = "forward")
	rda_reduced_b <- ordistep(rda_reduced, scope = formula(rda_0), direction = "backward")

	# Build a model manually
	rda_m <- update(rda_0, . ~ . + Mean.Temp + slope_gradient + Northness + Eastness) 
	head(summary(rda_m))
	anova(rda_m, by = "terms")

}

