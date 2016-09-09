library(vegetarian)
library(data.table)

setwd("H:/My Documents/PhD Research PFR folder/LBI_miseq_analysis_stuff/LBI_U8_OTUs_analyses/")
#setwd("Y:/PhD/PhD Research/NZGL01401_analysis/")

### Load OTU datasets ###
get.datasets <- function(genes){
  df.list <- list()
  n <- 1
  for(g in genes){
    #print(g)
    #f <- Sys.glob(paste0("OTUtables/", g, "*otutable.txt"))
    f <- Sys.glob(paste0("OTUtables/", g, "*otutable_by_plot.txt"))
    df <- read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    if(g == "ShCO1"){ # Fix COI names
      g <- "COI-300"
    }
    if(g == "FolCO1"){
      g <- "COI-650"
    }
    df.list[[n]] <- df
    label <- paste(g)
    print(label)
    print(nrow(df))
    names(df.list)[[n]] <- label
    n <- n + 1
  }
  return(df.list)
}

### Get overall alpha diversity stats per gene ###
get.overall.div <- function(df.list){
  all.div <- list()
  n <- 1
  for(df in df.list){
    div.list.1 <- list()
    div.list.2 <- list()
    print(names(df.list)[[n]])
    df.1 <- t(df)
    ### Min 1 stats ###
    a0 <- d(df.1, q = 0, lev = "alpha")
    a1 <- d(df.1, q = 1, lev = "alpha")
    a2 <- d(df.1, q = 2, lev = "alpha")
    b0 <- d(df.1, q = 0, lev = "beta")
    b1 <- d(df.1, q = 1, lev = "beta")
    b2 <- d(df.1, q = 2, lev = "beta")
    g0 <- d(df.1, q = 0, lev = "gamma")
    g1 <- d(df.1, q = 1, lev = "gamma")
    g2 <- d(df.1, q = 2, lev = "gamma")
    div.1 <- list("min" = "min1", "gene" = names(df.list)[[n]], 
                  "alpha0" = a0, "alpha1" = a1, "alpha2" = a2,
                  "beta0" = b0, "beta1" = b1, "beta2" = b2,
                  "gamma0" = g0, "gamma1" = g1, "gamma2" = g2)
    div.list.1[[n]] <- div.1
    
    ### Min 2 stats ###
    df.2 <- t(df.1[rowSums(df.1[,1:28]) > 1, 1:28])
    a0 <- d(df.2, q = 0, lev = "alpha")
    a1 <- d(df.2, q = 1, lev = "alpha")
    a2 <- d(df.2, q = 2, lev = "alpha")
    b0 <- d(df.2, q = 0, lev = "beta")
    b1 <- d(df.2, q = 1, lev = "beta")
    b2 <- d(df.2, q = 2, lev = "beta")
    g0 <- d(df.2, q = 0, lev = "gamma")
    g1 <- d(df.2, q = 1, lev = "gamma")
    g2 <- d(df.2, q = 2, lev = "gamma")
    div.2 <- list("min" = "min2", "gene" = names(df.list)[[n]], 
                  "alpha0" = a0, "alpha1" = a1, "alpha2" = a2,
                  "beta0" = b0, "beta1" = b1, "beta2" = b2,
                  "gamma0" = g0, "gamma1" = g1, "gamma2" = g2)
    div.list.2[[n]] <- div.2
    all.div[[n]] <- rbind(rbindlist(div.list.1), rbindlist(div.list.2))
    n <- n + 1
  }
  return(rbindlist(all.div))
}  

###############################################################################
genes <- c("16S","18S", "26S", "ITS", "FolCO1", "ShCO1")
df.list <- get.datasets(genes)
overall.div <- get.overall.div(df.list)

write.table(overall.div, file = "Overall_div_stats_by_gene_by_plot.txt", sep = "\t", quote = FALSE, col.names = NA)
