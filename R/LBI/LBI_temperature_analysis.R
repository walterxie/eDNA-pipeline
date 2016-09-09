library(ggplot2)
library(data.table)
library(reshape2)

setwd("D:/PhD_folder/LBI_miseq_analysis_stuff/LBI_environmental_data/")

df <- read.table("LBI_all_temperature_data.txt", sep="\t", header = T, check.names = F)

envdata_byplot <- read.table("LBI_all_env_data_by_plot.txt", sep="\t", header = T, check.names = F)
envdata_bysubplot <- read.table("LBI_all_env_data_by_subplot.txt", sep="\t", header = T, check.names = F)

# Fix sample names
colnames(df) <- gsub("Plot 9L", "CM30c30 L", colnames(df))
colnames(df) <- gsub("Plot 10E", "LB1 E", colnames(df))
colnames(df) <- gsub("\\s\\(.*", "", colnames(df)) # remove names in brackets
colnames(df) <- gsub("\\s([0-9])([A-Z])", "0\\1-\\2", colnames(df))
colnames(df) <- gsub("\\s([A-Z])", "-\\1", colnames(df))

# Separate surface level and 10cm depth measurements
d0cm <- df[, grep("\\s0cm", colnames(df))] 
d0cm$date.time <- df$"plot/subplot/depth"
d10cm <- df[, grep("-10cm", colnames(df))] 
d10cm$date.time <- df$"plot/subplot/depth"

colnames(d0cm) <- gsub("-[A-Z]\\s0cm", "", colnames(d0cm))
colnames(d10cm) <- gsub("-[A-Z]\\s-10cm", "", colnames(d10cm))

### Calculate overall annual mean temperature for each site
d0.ann.mean <- as.data.frame(colMeans(d0cm[, 1:28], na.rm = TRUE))
d10.ann.mean <- as.data.frame(colMeans(d10cm[, 1:28], na.rm = TRUE))
colnames(d0.ann.mean)[1] <- "mean"
colnames(d10.ann.mean)[1] <- "mean"

# Match plots with elevation
d0.ann.mean$Elevation <- envdata_byplot$Elevation[match(rownames(d0.ann.mean), envdata_byplot$Plot)]
d10.ann.mean$Elevation <- envdata_byplot$Elevation[match(rownames(d10.ann.mean), envdata_byplot$Plot)]

# Plot annual mean temperature vs. elevation
par(mfrow=c(2,2))
plot(d0.ann.mean$mean ~ d0.ann.mean$Elevation)
plot(d10.ann.mean$mean ~ d10.ann.mean$Elevation)

### Calculate annual mean temperature at time x for each site
d0.ann.mean.12pm <- as.data.frame(colMeans(d0cm[grepl("12:00:00 PM", d0cm$date.time), 1:28], na.rm = TRUE))
d10.ann.mean.12pm <- as.data.frame(colMeans(d10cm[grepl("12:00:00 PM", d0cm$date.time), 1:28], na.rm = TRUE))

colnames(d0.ann.mean.12pm)[1] <- "mean"
colnames(d10.ann.mean.12pm)[1] <- "mean"

# Match plots with elevation
d0.ann.mean.12pm$Elevation <- envdata$Elevation[match(rownames(d0.ann.mean.12pm), envdata$Plot)]
d10.ann.mean.12pm$Elevation <- envdata$Elevation[match(rownames(d10.ann.mean.12pm), envdata$Plot)]

# Plot annual mean temperature at time x vs. elevation
par(mfrow=c(1,2))
plot(d0.ann.mean.12pm$mean ~ d0.ann.mean.12pm$Elevation)
plot(d10.ann.mean.12pm$mean ~ d10.ann.mean.12pm$Elevation)

### Add means to envdata
envdata_byplot$mean_T_surface <- d0.ann.mean$mean[match(envdata_byplot$Plot, rownames(d0.ann.mean))]
envdata_byplot$mean_T_10cm <- d10.ann.mean$mean[match(envdata_byplot$Plot, rownames(d10.ann.mean))]
envdata_bysubplot$mean_T_surface <- d0.ann.mean$mean[match(envdata_bysubplot$Plot, rownames(d0.ann.mean))]
envdata_bysubplot$mean_T_10cm <- d10.ann.mean$mean[match(envdata_bysubplot$Plot, rownames(d10.ann.mean))]

envdata_byplot$mean_T_surface_12pm <- d0.ann.mean.12pm$mean[match(envdata_byplot$Plot, rownames(d0.ann.mean.12pm))]
envdata_byplot$mean_T_10cm_12pm <- d10.ann.mean.12pm$mean[match(envdata_byplot$Plot, rownames(d10.ann.mean.12pm))]
envdata_bysubplot$mean_T_surface_12pm <- d0.ann.mean.12pm$mean[match(envdata_bysubplot$Plot, rownames(d0.ann.mean.12pm))]
envdata_bysubplot$mean_T_10cm_12pm <- d10.ann.mean.12pm$mean[match(envdata_bysubplot$Plot, rownames(d10.ann.mean.12pm))]

plot(envdata_byplot$mean_T_10cm ~ envdata_byplot$mean_T_surface)
plot(envdata_byplot$mean_T_10cm_12pm ~ envdata_byplot$mean_T_surface_12pm)

write.table(envdata_byplot, file = "LBI_all_env_data_by_plot.txt", sep="\t", col.names = NA, row.names = TRUE, quote = F)
write.table(envdata_bysubplot, file = "LBI_all_env_data_by_subplot.txt", sep="\t", col.names = NA, row.names = TRUE, quote = F)

### Plot all measurements together
## Surface measurements
d0.l <- melt(d0cm, id.vars = c("date.time"))
d0.l$date <- sapply(strsplit(as.character(d0.l$date.time), "\\s"), "[[", 1)
d0.l$time <- paste(sapply(strsplit(as.character(d0.l$date.time), "\\s"), "[[", 2),
                      sapply(strsplit(as.character(d0.l$date.time), "\\s"), "[[", 3))
d0.l$month <- sapply(strsplit(as.character(d0.l$date), "/"), "[[", 1)
d0.l$Elevation <- envdata$Elevation[match(d0.l$variable, envdata$Plot)] 

ggplot(na.omit(d0.l)) + geom_jitter(aes(x = factor(Elevation), y = value, group=month, colour = month)) 
+ 
  facet_grid(~month)
ggplot(na.omit(d0.l)) + geom_boxplot(aes(x = factor(Elevation), y = value))
ggplot(na.omit(d0.l)) + geom_boxplot(aes(x = factor(Elevation), y = value)) + facet_wrap(~month)

# Plot all measurements from time x together vs. elevation
d0.sub <- na.omit(d0.l[grepl("12:00:00 PM", d0.l$time), ])
summary(d0.sub)

ggplot(na.omit(d0.sub)) + geom_jitter(aes(x = factor(Elevation), y = value, group=month, colour = month)) + 
      facet_grid(~month)
ggplot(na.omit(d0.sub)) + geom_smooth(aes(x = factor(Elevation), y = value, group=month)) + 
      facet_wrap(~month)

### 10cm depth measurements
d10.l <- melt(d10cm, id.vars = c("date.time"))
d10.l$date <- sapply(strsplit(as.character(d10.l$date.time), "\\s"), "[[", 1)
d10.l$time <- paste(sapply(strsplit(as.character(d10.l$date.time), "\\s"), "[[", 2),
                   sapply(strsplit(as.character(d10.l$date.time), "\\s"), "[[", 3))
d10.l$month <- sapply(strsplit(as.character(d10.l$date), "/"), "[[", 1)
d10.l$Elevation <- envdata$Elevation[match(d10.l$variable, envdata$Plot)] 

ggplot(na.omit(d10.l)) + geom_point(aes(x = factor(Elevation), y = value, group=month, colour = Elevation)) + 
  facet_grid(~month)
ggplot(na.omit(d10.l)) + geom_boxplot(aes(x = factor(Elevation), y = value))
ggplot(na.omit(d10.l)) + geom_boxplot(aes(x = factor(Elevation), y = value)) + facet_wrap(~month)
ggplot(na.omit(d10.l)) + geom_smooth(aes(x = factor(Elevation), y = value, group=month), method="lm") + 
  facet_wrap(~month)

# Plot all measurements from time x together vs. elevation
d10.sub <- na.omit(d10.l[grepl("12:00:00 PM", d10.l$time), ])
summary(d10.sub)

ggplot(na.omit(d10.sub)) + geom_jitter(aes(x = factor(Elevation), y = value, group=month, colour = month)) + 
  facet_grid(~month)
ggplot(na.omit(d10.sub)) + geom_smooth(aes(x = factor(Elevation), y = value, group=month)) + 
  facet_wrap(~month)


### Use data.table to get mean by month?
dt <- data.table(d10.l)

dt[, mean, by=month] 

#################################
#library(lubridate)
#ydm_hms(d.0cm.l$date.time, tz = "Pacific/Auckland")

model = lm(d0.sub$value ~ d0.sub$Elevation)
summary(model)
par(mfrow=c(2,2))
plot(model)

plot(d0.sub$Elevation, d0.sub$value, xlim=c(0,700), ylim=c(0,30))
abline(model)

model = lm(d0.sub$value ~ d0.sub$Elevation + I(d0.sub$Elevation^3) + d0.sub$date)
summary(model)
par(mfrow=c(2,2))
plot(model)

plot(loess(d0.sub$value ~ d0.sub$Elevation))
abline(model)

model = lm(d0.sub$value ~ d0.sub$Elevation + d0.sub$date)
summary(model)
par(mfrow=c(2,2))
plot(model)

y <- d0.sub$value
x <- d0.sub$Elevation
x2 <- d0.sub$date



