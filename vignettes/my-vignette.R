## ------------------------------------------------------------------------
library(crossmeta)
library(ccmap)
library(lydata)
data_dir <- system.file("extdata", package = "lydata")

# gather all GSEs
gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")

# load previous crossmeta differential expression analysis
anals <- load_diff(gse_names, data_dir)

# run meta-analysis
es <- es_meta(anals)

# extract moderated adjusted standardized effect sizes
dprimes <- get_dprimes(es)

# query signature
query_sig <- dprimes$meta

## ------------------------------------------------------------------------
library(ccdata)

# load drug signatures
data(cmap_es)

## ------------------------------------------------------------------------
library(ccmap)

# query drug signatures using all common genes
top_drugs <- query_drugs(query_sig)

# query drug signatures using a range of query gene sizes
# range_res <- range_query_drugs(query_sig)

head(top_drugs[[1]])


## ----eval=FALSE----------------------------------------------------------
#  library(ccmap)
#  
#  # make drug combination data base
#  make_drug_combos()
#  
#  db_dir <- "/path/to/drug_combos.sqlite"
#  
#  # query drug combination signatures using all common genes
#  top_combos <- query_combos(query_sig, db_dir)
#  
#  # query drug signatures using a range of query gene sizes
#  # (takes approximately 5 hours using a Intel® Core™ i7-6700K)
#  range_res <- range_query_combos(query_sig, db_dir)
#  

