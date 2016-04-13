library(RSQLite)
library(dplyr)

source("~/Documents/Batcave/GEO/1-meta/MAMA_utils.R")
source("~/Documents/Batcave/GEO/2-cmap/cmap_utils.R")


#---------
# SETUP
#---------

setwd("/home/alex/Documents/Batcave/GEO/2-cmap/data/processed/es")

es2 <- readRDS("~/Documents/Batcave/GEO/1-meta/data/MAMA/cr/es2.rds")
drug_es <- readRDS("~/Documents/Batcave/GEO/2-cmap/data/processed/es/genes_drug_es.rds")

query_genes <- get_zscores(es2)$meta
drug_genes_n <- 500

#get up/dn query genes
query_genes <- setup_query_genes(query_genes, drug_info=drug_es)

#ordered combo names allow fast selection from database
names_combos <- readRDS("~/Documents/Batcave/GEO/2-cmap/data/processed/es/names_combos.rds")
names_probes <- readRDS("~/Documents/Batcave/GEO/2-cmap/data/processed/es/names_probes.rds")
ncombos <- length(names_combos)

#get annotation info
library("hgu133a.db")
map <- AnnotationDbi::select(hgu133a.db, names_probes, "SYMBOL")
map <- map[!is.na(map$SYMBOL),]

#df for top drugs
cols <- c("overlap", "cross", "net", "drug_genes_n", "net_drug_pct")
top_combos <- data.frame(matrix(nrow=ncombos, 
                                ncol=5, dimnames=list(names_combos, cols)))

top_combos$drug_genes_n <- drug_genes_n
rm(es2, drug_es, names_combos, names_probes, cols)
gc()

#connect to db
db <- dbConnect(SQLite(), dbname="probes_drug_combos.sqlite")


#query db
pb <- txtProgressBar(min=492, max=654, style=3)

#---------
# QUERY
#---------

for (i in 492:654) {
  #get hypothesized effect size for drug combo
  statement <- paste("SELECT * from combo_tstats WHERE rowid BETWEEN", (i*1309)-1308, "AND", i*1309)
  combos_es <- dbGetQuery(db, statement)

  #format result
  drug_combos <- combos_es$drug_combo

  row.names(combos_es) <- drug_combos
  combos_es <- as.data.frame(t(combos_es[,-1]))
  
  #annotate result
  combos_es <- combos_es[map$PROBEID, ]
  combos_es$SYMBOL <- map$SYMBOL

  for (combo in drug_combos) {

    combo_es <- combos_es[, c(combo, "SYMBOL")]
    colnames(combo_es) <- c("tstat", "SYMBOL")

    combo_es %>%
      group_by(SYMBOL) %>%
      arrange(desc(abs(tstat))) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      arrange(desc(tstat)) ->
      combo_es

    class(combo_es) <- "data.frame"
    row.names(combo_es) <- combo_es$SYMBOL
    combo_es <- combo_es[,-2, drop=F]

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