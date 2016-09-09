library(ggplot2)
library(data.table)
#source("H:/My Documents/PhD Research PFR folder/R stuff/R_LBI_Miseq/IWantHue.R")
theme_set(theme_bw(base_size=8))

###############################################################################
# Compare composition/structure of old/new communities (requires OTU identifications to be retained between networks)
setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/Network_analysis/")

#outfile <- "Network_analysis_20-80_community_consistency.txt"

comms_a <- read.table("Networks_by_OTUs_all_genes_veg_v3_phy_20-80_b/LBI_all_genes_veg_ggnet_by_OTUs_fg_all_communities_vertices_data.txt", sep = "\t", header = T)
comms_b <- read.table("Networks_by_OTUs_all_genes_veg_v3_phy_20-80_d/LBI_all_genes_veg_ggnet_by_OTUs_fg_all_communities_vertices_data.txt", sep = "\t", header = T)

#cat("\n\n\nNetwork community comparisons, 20-80 a. vs b.\n\n", file = outfile, append = TRUE)
#cat("\n\n\nNetwork community comparisons, 20-80 b. vs d.\n\n", file = outfile, append = TRUE)

### Check how many communities in each dataset ###
max(comms_a$comm)
max(comms_b$comm)

### Calculate actual results ###
### Set up matrix to hold data ###
m <- matrix(nrow = max(comms_a$comm), ncol = max(comms_b$comm))
  
r_names <- list() # to hold rownames
c_names <- list() # to hold column names
  
for(i in 1:max(comms_a$comm)){
  comm_a <- comms_a[comms_a$comm == i, ]
  nodes_a <- unique(comm_a$X)
  r_names[i] <- paste0("comm.", i, ".", length(nodes_a))
  for (j in 1:max(comms_b$comm)) {
    comm_b <- comms_b[comms_b$comm == j, ]
    nodes_b <- unique(comm_b$X)
    c_names[j] <- paste0("comm.", j, ".", length(nodes_b))
    #rownames(m)$j <- paste(length(nodes_b)) 
    matches <- match(nodes_a, nodes_b)
    m[i,j] <- length(matches[!is.na(matches)])
    print(paste(i, "nodes_a:", length(nodes_a), j , "nodes_b:", length(nodes_b), "matching:", length(matches[!is.na(matches)])))
  }
}
colnames(m) <- c_names
rownames(m) <- r_names
actual <- melt(m)

### Generate randomised results ###
results <- list()
for(n in 1:1000){
  print(paste("randomising", n))
  comms_a$comm <- sample(comms_a$comm)
  comms_b$comm <- sample(comms_b$comm)

  ### Set up matrix to hold data ###
  m <- matrix(nrow = max(comms_a$comm), ncol = max(comms_b$comm))
  
  r_names <- list() # to hold rownames
  c_names <- list() # to hold column names
  
  for(i in 1:max(comms_a$comm)){
    comm_a <- comms_a[comms_a$comm == i, ]
    nodes_a <- unique(comm_a$X)
    r_names[i] <- paste0("comm.", i, ".", length(nodes_a))
    for (j in 1:max(comms_b$comm)) {
      comm_b <- comms_b[comms_b$comm == j, ]
      nodes_b <- unique(comm_b$X)
      c_names[j] <- paste0("comm.", j, ".", length(nodes_b))
      #rownames(m)$j <- paste(length(nodes_b)) 
      matches <- match(nodes_a, nodes_b)
      m[i,j] <- length(matches[!is.na(matches)])
      #print(paste(i, "nodes_a:", length(nodes_a), j , "nodes_b:", length(nodes_b), "matching:", length(matches[!is.na(matches)])))
    }
  }
  colnames(m) <- c_names
  rownames(m) <- r_names
  
  results[[n]] <- melt(m)
}

res <- rbindlist(results)

res$type <- "predicted"
actual$type <- "actual"
all.data <- rbind(res, actual)
all.data$combo <- paste0(all.data$Var1, all.data$Var2)

get.color <- function(z){
  d <- all.data[grep(z, all.data$combo), ]
  if(d$value[d$type == "actual"] > max(d$value[d$type == "predicted"])){
    #if(d$type[grep(max(d$value), d$value)] == "actual"){
    d$col <- "Actual"
  }else if(d$value[d$type == "actual"] == max(d$value[d$type == "predicted"])){
    d$col <- "Same"
  }else{
    d$col <- "Predicted"
  }
  return(d)
}  

new.data <- data.frame()
for(z in unique(all.data$combo)){
  d <- get.color(z)
  new.data <- rbind(new.data, d)
}

ggplot(new.data, aes(x = value, fill = col)) + geom_histogram(binwidth = 1) +
  #ggplot(res, aes(x = value)) + geom_freqpoly(binwidth = 1) +
  geom_vline(data = actual, aes(xintercept = value)) + 
  scale_fill_manual(values = c("#e41a1c","grey","#377eb8")) +
  facet_grid(Var1 ~ Var2, scales = "free") +
  theme(panel.grid = element_blank(), legend.position="none")

zzz <- all.data[all.data$Var2 == "comm.3.316",]
write.table(zzz, file = "blahblah.txt", sep = "\t", quote = F, col.names = NA)
