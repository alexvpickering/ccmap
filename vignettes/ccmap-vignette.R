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
top_drugs <- query_drugs(query_sig)

# correctly identifies LY-294002 as best match among drug signatures
# other PI3K inhibitors are also identified among top matching drugs
head(top_drugs)

## ---- message=FALSE, warning=FALSE---------------------------------------
# query all 856086 combinations (takes ~10 minutes on Intel Core i7-6700)
# top_combos <- query_combos(query_sig)

# query only combinations with LY-294002
top_combos <- query_combos(query_sig, include='LY-294002', ncores=1)

## ---- message=FALSE, warning=FALSE---------------------------------------
# Times on Intel Core i7-6700 with MRO+MKL
# requires ~8-10GB of RAM

method  <- 'ml'
include <- names(head(top_drugs))

# query all 856086 combinations (~2 hours)
# top_combos <- query_combos(query_sig, method)

# query combinations with top single drugs (~1 minute)
# top_combos <- query_combos(query_sig, method, include)

sessionInfo()

