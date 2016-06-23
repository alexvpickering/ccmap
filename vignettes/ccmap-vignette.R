## ---- message=FALSE, warning=FALSE---------------------------------------
library(crossmeta)
library(ccmap)

# microarray data from studies using drug LY-294002
library(lydata)
data_dir <- system.file("extdata", package = "lydata")

# gather all GSEs
gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")

# load previous crossmeta differential expression analysis
anals <- load_diff(gse_names, data_dir)

# run meta-analysis
es <- es_meta(anals)

# contribute your signature to our public meta-analysis database
# contribute(anals, subject = "LY-294002")

# extract moderated adjusted standardized effect sizes
dprimes <- get_dprimes(es)

# query signature
query_sig <- dprimes$meta

## ---- message=FALSE, warning=FALSE---------------------------------------
library(ccdata)

# load drug signatures
data(cmap_es)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(ccmap)
top_drugs <- query_drugs(query_sig)

# correctly identifies LY-294002 as best match among drug signatures
# other PI3K inhibitors are also identified among top matching drugs
head(top_drugs)

## ----eval=FALSE----------------------------------------------------------
#  library(ccmap)
#  # this takes ~6 minutes per drug in 'with'
#  drug_combos <- predict_combos(with = c("LY-294002", "dilazep"))
#  
#  # query drug combination signatures using a range of query gene sizes
#  range_res <- query_drugs(query_sig, drug_combos)

