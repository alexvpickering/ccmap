library(ggplot2)
library(reshape2)

source("~/Documents/Batcave/GEO/2-cmap/cytoscape_utils.R")
source("~/Documents/Batcave/GEO/2-cmap/eval_atc.R")

#---------------

setwd("~/Documents/Batcave/GEO/2-cmap/data/processed")

#load in cmap data
drug_es <- readRDS("~/Documents/Batcave/GEO/2-cmap/data/processed/es/genes_drug_es.rds")
atc4 <- readRDS("~/Documents/Batcave/GEO/2-cmap/data/atc/atc4.rds")
drugs <- colnames(drug_es)

#get similarity matrix long
tsim <- get_tan_sim(drug_es, es=T)

#get tp/fp rates
trates <- get_rates(tsim, atc4)

#put tp/fp rates into dataframe
trates_df <- data.frame(fp=trates$fp, tp=trates$tp)

#plot tp vs fp
ggplot(trates_df, aes(fp, tp)) + geom_point()