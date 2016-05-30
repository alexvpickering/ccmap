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
# range_res <- query_drugs(query_sig, step = 100)

head(top_drugs)

## ----eval=FALSE----------------------------------------------------------
#  library(ccmap)
#  # this takes ~12 minutes per drug in 'with'
#  drug_combos <- predict_combos(with = c("LY-294002", "dilazep"))
#  
#  # query drug combination signatures using a range of query gene sizes
#  range_res <- query_drugs(query_sig, drug_combos, step = 100)

