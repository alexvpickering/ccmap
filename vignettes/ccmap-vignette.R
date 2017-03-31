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
# load drug signatures
data(l1000_es)

## ---- message=FALSE, warning=FALSE---------------------------------------
top_cmap  <- query_drugs(query_sig, cmap_es)
top_l1000 <- query_drugs(query_sig, l1000_es)

# LY-294002 best match among 1309 cmap signatures
# other PI3K inhibitors are also identified among top matching drugs
head(top_cmap, 4)

# LY-294002 matches 4 of top 10 l1000 signatures (230,829 total) 
# other PI3K inhibitors are also identified among top matching drugs
head(top_l1000, 4)

## ---- message=FALSE, warning=FALSE---------------------------------------
# remove genes in cmap_es that are not measured in l1000_es
cmap_lm <- cmap_es[row.names(l1000_es), ]

# query using genes common to cmap_es and l1000_es
top_cmap_lm  <- query_drugs(query_sig, cmap_lm)

## ---- message=FALSE, warning=FALSE---------------------------------------
# query all 856086 combinations (takes ~2 minutes on Intel Core i7-6700)
# top_combos <- query_combos(query_sig, cmap_es)

# query only combinations with LY-294002
top_combos <- query_combos(query_sig, cmap_es, include='LY-294002', ncores=1)

## ---- message=FALSE, warning=FALSE---------------------------------------

# query only combinations with LY-294002_NKDBA_10um_24h
top_combos <- query_combos(query_sig, l1000_es, include='LY-294002_NKDBA_10um_24h', ncores=1)

# query combinations with all LY-294002 signatures
# top_combos <- query_combos(query_sig, l1000_es, include='LY-294002')

## ---- message=FALSE, warning=FALSE---------------------------------------
# Times on Intel Core i7-6700 with MRO+MKL
# requires ~8-10GB of RAM

method  <- 'ml'
include <- names(head(top_cmap))

# query all 856086 combinations (~2 hours)
# top_combos <- query_combos(query_sig, 'cmap', method)

# query combinations with top single drugs (~1 minute)
# top_combos <- query_combos(query_sig, 'cmap', method, include)

sessionInfo()

