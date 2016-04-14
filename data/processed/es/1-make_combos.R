library(RSQLite)
library(randomForestSRC)

setwd("/home/alex/Documents/Batcave/GEO/2-cmap/data/processed/es")

#---------

#-------
# SETUP
#-------


#load model & cmap data
mod1        <- readRDS("~/Documents/Batcave/GEO/2-cmap/data/combos/mod1.rds")
mod2        <- readRDS("~/Documents/Batcave/GEO/2-cmap/data/combos/mod2.rds")
cmap_probes <- readRDS("~/Documents/Batcave/GEO/2-cmap/data/processed/es/probes_top_tables.rds")
cmap_array  <- readRDS("~/Documents/Batcave/GEO/2-cmap/data/combos/cmap_array.rds")
probes      <- mod1$probes
drugs       <- cmap_probes$drugs


#probes_data <- cmap_probes$data[probes, ]

#indices for probes_data
drug_inds <- sapply(seq_along(drugs), function(x) c(x*7-6, x*7))
colnames(drug_inds) <- drugs

#make blank table in SQLite database
db.pdc <- dbConnect(SQLite(), dbname="genes_drug_combos.sqlite")
cols <- paste(shQuote(probes), "FLOAT", collapse=", ")

statement <- paste("CREATE TABLE combo_preds",
                   "(drug_combo TEXT, ", cols, ")", sep="")

dbSendQuery(db.pdc, statement)


#get list of unique drug combos
combo_pairs <- combn(drugs, 2, simplify=F)
combo_names <- sapply(combo_pairs, function(x) paste(x[1], x[2], sep=" + "))
names(combo_pairs) <- combo_names

#store names for fast retrieval from database
saveRDS(combo_names, "combo_names.rds")

#-----------------------
# POPULATE combo_preds
#-----------------------
# 856086 combos = 1309 * 654

combine_arrays <- function(combos, cmap_array) {
	
	drugs1 <- sapply(combos, function(x) x[[1]])
	drugs2 <- sapply(combos, function(x) x[[2]])

  data <- cmap_array$data
  cols <- cmap_array$col_names

	#bind drug1 and drug2 data by columns
	drugs12 <- list()
	for (i in seq_along(combos)) {

    drug1_data <- data[[drugs1[i]]]
    drug2_data <- data[[drugs2[i]]]

		drugs12[[i]] <- c(drug1_data, drug2_data)
	}

	#bind drug12s data by rows
  drugs12 <- t(as.data.frame(drugs12))

  row.names(drugs12) <- names(combos)
  colnames(drugs12)  <- cols

  return(as.data.frame(drugs12))
}




#~10hr and 180Gb
for (i in 1:1309) {

	#get array data for 654 combos
  combos <- combo_pairs[(i*654-653):(i*654)]
	combos_data <- combine_arrays(combos, cmap_array)

	#make predictions with mod1
  start <- Sys.time()
	preds1 <- predict.rfsrc(mod1$model, combos_data)
  end <- Sys.time()
  end-start

	#get probe data for 654 combos

	#add predictions from mod1

	#make predictions with mod2


    combo_tstat <- es_probes[, d1] + es_probes[, d2]
    combo_tstat <- paste(combo_tstat, collapse=", ")

    statement <- paste("INSERT INTO combo_tstats VALUES (", shQuote(combo_name), sep="")
    statement <- paste(statement, ", ", combo_tstat, ")", sep="")

    dbSendQuery(conn = db.pdc, statement)
}





  #--------------
  # PROBE COMBOS
  #--------------

  #df for predicted probabilities of combos
  probe_preds <- data.frame(row.names=probes)

  for (drug2 in other_drugs) {

    #get cmap probe data for drug2
    drug2_data <- cmap_data[, grepl(drug2, colnames(cmap_data))]
    colnames(drug2_data) <- gsub(drug2, "drug2", colnames(drug2_data))

    #get mod1 preds
    combo_data1 <- by_array(top_data[mod1$order, ], drug2_data[mod1$order, ])
    preds1 <- predict.rfsrc(mod1$model, combo_data1)



    #bind top and drug2 then get preds for combo
    combo_data <- cbind(top_data, drug2_data)
    probe_preds[, drug2] <- predict(model, as.matrix(combo_data))
  }

  #-----------------
  # PROBES to GENES
  #-----------------

  #annotate probe preds with SYMBOL
  suppressMessages(library("hgu133a.db"))
  suppressMessages(map <- AnnotationDbi::select(hgu133a.db, probes, "SYMBOL"))
  map <- map[!is.na(map$SYMBOL),]

  probe_preds <- probe_preds[map$PROBEID, ]
  probe_preds$SYMBOL <- map$SYMBOL

  #df for predicted probabilities of combos at gene level
  gene_preds <- data.frame(row.names=unique(map$SYMBOL))

  for (drug2 in other_drugs) {

    #for duplicated genes, choose probability furthest from 0.5
    combo_preds <- probe_preds[, c(drug2, "SYMBOL")]
    colnames(combo_preds) <- c("pred", "SYMBOL")
    combo_preds$dist <- combo_preds$pred - 0.5

    combo_preds %>%
      group_by(SYMBOL) %>%
      arrange(desc(abs(dist))) %>%
      dplyr::slice(1) %>%
      ungroup() ->
      combo_preds

    #add to gene_preds
    class(combo_preds) <- "data.frame"
    gene_preds[combo_preds$SYMBOL, drug2] <- combo_preds$dist
  }

  #------------
  # BEST COMBO
  #------------

  #get best combo of top_drugs and other_drugs
  best_combo <- get_top_drugs(query_genes, query_n, as.matrix(gene_preds),
                              drug_genes_n=drug_genes_n, es=T)$table[1, ]

  combo_table <- rbind(combo_table, best_combo)
  return(combo_table)
}







#load in matrix of probe-level t-statistics (log2FC/SD)
es_probes <- readRDS("probes_drug_es.rds")

#make blank table in SQLite database
db.pdc <- dbConnect(SQLite(), dbname="probes_drug_combos.sqlite")

probes <- row.names(es_probes)
drugs <- colnames(es_probes)
cols <- paste(shQuote(probes), "FLOAT", collapse=", ")

statement <- paste("CREATE TABLE combo_tstats",
                   "(drug_combo TEXT, ", cols, ")", sep="")

dbSendQuery(db.pdc, statement)



#get list of unique drug combos
combo_pairs <- combn(drugs,2, simplify=F)
combo_names <- sapply(combo_pairs, function(x) paste(x[1], x[2], sep=" + "))
names(combo_pairs) <- combo_names

#store names for fast retrieval from database
saveRDS(combo_names, "combo_names.rds")


#-----------------------
# POPULATE combo_tstats
#-----------------------

#~10hr and 180Gb
for (i in seq_along(combo_pairs)) {
    combo_name <- names(combo_pairs[i])
    
    d1 <- combo_pairs[[i]][1]
    d2 <- combo_pairs[[i]][2]

    combo_tstat <- es_probes[, d1] + es_probes[, d2]
    combo_tstat <- paste(combo_tstat, collapse=", ")

    statement <- paste("INSERT INTO combo_tstats VALUES (", shQuote(combo_name), sep="")
    statement <- paste(statement, ", ", combo_tstat, ")", sep="")

    dbSendQuery(conn = db.pdc, statement)
}



