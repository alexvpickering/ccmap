#----------
# APPROACH
#----------

# have t-statistic (log2FC/SD) - - > assume additive at FC level
#
#
#               | control   |   drug-1    |   drug-2   |   drug-1 + drug-2 (additive model)
#----------------------------------------------------------------------------------
# abs (FC)      |   10      |   20 (2x)    |   20 (2x)   |       40 (4x)    
# log2FC        |   N/A     |   1          |   1         |       2


#-------
# SETUP
#-------

library(RSQLite)

setwd("/home/alex/Documents/Batcave/GEO/2-cmap/data/processed/es")

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



