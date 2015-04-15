# validate input
if(!exists("inputEnvData")) stop("inputEnvData is missing !")

envData <- read.table(inputEnvData, header=T, row.names=1, sep="\t")

if (exists("adonis_formula")) {

	aov <- adonis(formula = eval(adonis_formula) ~ pH + Olsen.P + EC + Organic.C + 
				  Total.N + NO3.N + NH4.N + Water.Content, 
				  data = envData,
				  permutations = 1999, 
				  strata = as.factor(sampleFactor))

} else {
	aov <- adonis(d.hornMorisita~as.factor(sampleFactor), permutations=1999)
}


#LJ12027 <- read.table("../data/environmental/LJ12027.txt", 
#                      header=T, row.names=1, sep="\t")
#
#aov <- adonis(formula = d.brayBin ~ pH + Olsen.P + EC + Organic.C + 
#              Total.N + NO3.N + NH4.N + Water.Content, 
#              data = LJ12027,
#              permutations = 1999, 
#              strata = as.factor(sampleFactor))
#
#adonis(formula = d.turnover ~ LJ12027$pH + LJ12027$Olsen.P + LJ12027$EC + 
#       LJ12027$Organic.C + LJ12027$Total.N + LJ12027$NO3.N + LJ12027$NH4.N + 
#       LJ12027$Water.Content, 
#       permutations = 999, 
#       strata = as.factor(sampleFactor))
