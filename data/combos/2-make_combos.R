library(RSQLite)
library(data.table)
library(xgboost)
library(dplyr)

source("~/Documents/Batcave/GEO/2-cmap/combo_utils.R")
setwd("~/Documents/Batcave/GEO/2-cmap/data/combos")


#-------
# SETUP
#-------


#load model & cmap data
mod    <- readRDS("model.rds")
cmap   <- readRDS("~/Documents/Batcave/GEO/2-cmap/data/processed/es/probes_top_tables.rds")
data   <- cmap$data
drugs  <- cmap$drugs
probes <- row.names(data)

#indices for drugs
inds <- sapply(seq_along(drugs), function(x) c(x*7-6, x*7))
colnames(inds) <- drugs

#list of unique drug combos
combo_pairs <- combn(drugs, 2, simplify=F)
combo_names <- sapply(combo_pairs, function(x) paste(x[1], x[2], sep=" + "))
names(combo_pairs) <- combo_names

#load annotation information
suppressMessages(library("hgu133a.db"))
suppressMessages(map <- AnnotationDbi::select(hgu133a.db, probes, "SYMBOL"))
map <- map[!is.na(map$SYMBOL),]
genes <- unique(map$SYMBOL)

#make blank table in SQLite database
db.pdc <- dbConnect(SQLite(), dbname="genes_drug_combos.sqlite")
cols <- paste(shQuote(genes), "FLOAT", collapse=", ")

statement <- paste("CREATE TABLE combo_preds",
                                     "(drug_combo TEXT, ", cols, ")", sep="")

dbSendQuery(db.pdc, statement)



# 856086 combo_pairs = 1309 * 654
for (i in 1:1309) {

    #-------------------
    # PROBE PREDICTIONS
    #-------------------

    #get array data for 654 combos
    pairs <- combo_pairs[(i*654-653):(i*654)]
    combo_data <- combine_cmap_pairs(pairs, data, inds)

    #make predictions with mod
    probe_preds <- predict(mod, combo_data)
    
    #probe_preds vector to df
    dim(probe_preds) <- c(length(probes), length(pairs))
    colnames(probe_preds)  <- names(pairs)
    row.names(probe_preds) <- probes
    probe_preds <- as.data.frame(probe_preds)

    #-----------------
    # PROBES TO GENES
    #-----------------

    #annotate probe probe_preds with SYMBOL
    probe_preds <- probe_preds[map$PROBEID, ]
    probe_preds$SYMBOL <- map$SYMBOL

    probe_preds <- data.table(probe_preds)
    gene_preds <- probe_preds[, lapply(.SD, max), by=group]

    start=Sys.time()
    #df for predicted probabilities of combos at gene level
    gene_preds <- data.frame(row.names=unique(map$SYMBOL))

    for (combo in names(pairs)) {

        #for duplicated genes, choose probability furthest from 0.5
        combo_preds <- probe_preds[, c(combo, "SYMBOL")]
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
        gene_preds[combo_preds$SYMBOL, combo] <- combo_preds$dist
    }
    end=Sys.time()
    end-start

    #----------------
    # ADD TO SQL DB
    #----------------

    #format table
    gene_preds <- as.data.frame(t(gene_preds))
    gene_preds$drug_combo <- row.names(gene_preds)
    row.names(gene_preds) <- NULL

    dbWriteTable(db.pdc, "combo_preds", gene_preds, append = TRUE)

    statement <- paste("INSERT INTO combo_tstats VALUES (", shQuote(combo_name), sep="")
    statement <- paste(statement, ", ", combo_tstat, ")", sep="")

    dbSendQuery(conn = db.pdc, statement)
}



