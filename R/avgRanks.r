library(xtable)

# change config below
workingPath <- "~/Projects/NZGO/pilot2010/pipeline/"

inputT <- paste(workingPath, "plots-rank-table-gamma0.txt", sep="")
#inputT <- paste(workingPath, "plots-rank-table-beta1.txt", sep="")

rank_table <- read.table(inputT, header=T, row.names=1, check.names=FALSE)

rownames(rank_table) <- gsub("CM30C30", "Plot9", rownames(rank_table), ignore.case = T)
rownames(rank_table) <- gsub("LB1", "Plot10", rownames(rank_table), ignore.case = T)

rank_table <- rank_table[c(3:nrow(rank_table),1,2),]

print(xtable(rank_table))

# average table
avg_table <- data.frame(row.names=rownames(rank_table))

avg_table$eDNA <- rowMeans(rank_table[,1:6], na.rm=T)
avg_table$Trad <- rowMeans(rank_table[,7:10], na.rm=T)
avg_table$Trad.No.Birds <- rowMeans(rank_table[,7:9], na.rm=T)
avg_table$All <- rowMeans(rank_table, na.rm=T)
avg_table$All.No.Birds <- rowMeans(rank_table[,1:9], na.rm=T)

avg_table <- formatC(signif(as.matrix(avg_table),digits=2), digits=2,format="fg", flag="#")

# std table
std_table <- matrix(0,nrow=nrow(avg_table),ncol=ncol(avg_table))
colnames(std_table) <- colnames(avg_table)
rownames(std_table) <- rownames(avg_table)

for (r in 1:nrow(std_table)) {
	std_table[r,1] <- sd(rank_table[r,1:6], na.rm=T)
	std_table[r,2] <- sd(rank_table[r,7:10], na.rm=T)
	std_table[r,3] <- sd(rank_table[r,7:9], na.rm=T)
	std_table[r,4] <- sd(rank_table[r,], na.rm=T)
	std_table[r,5] <- sd(rank_table[r,1:9], na.rm=T)
}

std_table <- formatC(signif(std_table,digits=3), digits=3,format="fg", flag="#")

# final
final_table <- matrix(0,nrow=nrow(avg_table),ncol=ncol(avg_table))
colnames(final_table) <- colnames(avg_table)
rownames(final_table) <- rownames(avg_table)
for (r in 1:nrow(avg_table)) {
	final_table[r,] <- paste(avg_table[r,], " $\\pm$ ",  std_table[r,])
}

print(xtable(final_table),sanitize.text.function=function(x){x})

