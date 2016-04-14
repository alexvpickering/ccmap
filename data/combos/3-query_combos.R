library(RSQLite)
library(data.table)
library(foreach)
library(doMC)

source("~/Documents/Batcave/GEO/2-cmap/cmap_utils.R")
setwd("~/Documents/Batcave/GEO/2-cmap/data/combos")


#---------
# SETUP
#---------


#load CR dprimes
es2 <- readRDS("~/Documents/Batcave/GEO/1-meta/data/MAMA/cr/es2.rds")

#get query gene dprimes
query_genes <- get_dprimes(es2)$meta

#connect to db
db <- dbConnect(SQLite(), dbname="genes_drug_combos.sqlite")


#---------
# QUERY
#---------
# 856086 combos = 1309 * 109 * 6

registerDoMC(6)

res_list <- foreach(i=1:6) %dopar% {

  a <- i*109-108
  b <- i*109

  pb  <- txtProgressBar(min=a, max=b, style=3)
  res <- list()

  for (j in a:b) {

    #get preds for drug combos
    statement   <- paste("SELECT * from combo_tstats WHERE rowid BETWEEN", (j*1309)-1308, "AND", j*1309)
    combo_preds <- dbGetQuery(db, statement)

    #format result
    combo_names <- combo_preds$drug_combo
    combo_preds <- as.data.frame(t(combo_preds[,-1]))

    colnames(combo_preds)  <- combo_names

    #get top drug combos
    top_combos <- get_top_drugs(query_genes, drug_info=combo_preds, es=T)
    
    setTxtProgressBar(pb, j)
    res[[ length(res)+1 ]] <- top_combos
  }
  res <- rbindlist(res)
}






res <- rbindlist(res_list)


pb <- txtProgressBar(min=492, max=654, style=3)

for (i in 492:654) {



    #get top up/dn genes for combo
    combo_up <- head(row.names(combo_es), drug_genes_n/2)
    combo_dn <- tail(row.names(combo_es), drug_genes_n/2)

    #get overlap
    up <- intersect(combo_up, query_genes$up)
    dn <- intersect(combo_dn, query_genes$dn)
    tot <- length(up) + length(dn)

    #get cross overlap
    xup <- intersect(combo_up, query_genes$dn)
    xdn <- intersect(combo_dn, query_genes$up)
    xtot <- length(xup) + length(xdn)

    #add overlap/cross to top_combos
    top_combos[combo, "overlap"] <- tot
    top_combos[combo, "cross"] <- xtot
  }

  #cleanup
  saveRDS(top_combos, "top_combos4.rds")
  rm(combos_es, drug_combos, combo_es, combo_up,
     combo_dn, up, dn, tot, xup, xdn, xtot)

  gc()
  setTxtProgressBar(pb, i)
}