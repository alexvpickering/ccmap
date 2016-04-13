library(dplyr)

setwd("~/Documents/Batcave/GEO/2-cmap/data/processed/prl")
load("iorio_prls.ro")

#-------------------
#reference for order of probe names
probe_names <- iorio_prls[, 1]

#replace probe name with position of probe
for (i in seq_along(colnames(iorio_prls))) { 
  iorio_prls[,i] <- match(probe_names, iorio_prls[,i])
}

#set row names to probe names
row.names(iorio_prls) <- probe_names

#-----------
# ANNOTATE
#-----------

#map from probe names to SYMBOL
library(hgu133a.db)
map <- AnnotationDbi::select(hgu133a.db, probe_names, "SYMBOL")
map <- map[!is.na(map$SYMBOL), ]
iorio_prls <- iorio_prls[map$PROBEID, ] 

#where SYMBOL duplicated, keep SYMBOL with least central rank
drug_prls <- data.frame(SYMBOL=unique(map$SYMBOL), stringsAsFactors=F)
middle <- nrow(iorio_prls) / 2

drug_names <- names(iorio_prls)
iorio_prls[,"SYMBOL"] <- map$SYMBOL

for (drug in drug_names) {
  #get drug prl, SYMBOL, and distance from center
  drug_prl <- iorio_prls[, c(drug, "SYMBOL")]
  drug_prl$dist <- abs(drug_prl[, drug] - middle)
  
  drug_prl %>%
    group_by(SYMBOL) %>%
    arrange(desc(dist)) %>%
    slice(1) %>%
    ungroup() %>%
    select(-dist) %>%
    inner_join(drug_prls, by="SYMBOL") ->
    drug_prls
}

class(drug_prls) <- "data.frame"

for (drug in drug_names) {
  #re-rank (to fill in gaps)
  rank <- order(drug_prls[, drug])
  drug_prls[rank, drug] <- 1:nrow(drug_prls)
}

row.names(drug_prls) <- drug_prls[, "SYMBOL"]
drug_prls <- drug_prls[, drug_names]
drug_prls <- as.matrix(drug_prls)

saveRDS(drug_prls, file = "genes_drug_prls.rds")