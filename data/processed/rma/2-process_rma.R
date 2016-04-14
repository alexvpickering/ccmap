library(sva)
library(Biobase)
library(metaMA)

source("~/Documents/Batcave/GEO/1-meta/diff_utils.R")
setwd("~/Documents/Batcave/GEO/2-cmap/data/processed/rma")

#----------------
# LOAD IN DATA
#----------------

#load RMA processed data for each platform
ht_hga_ea <- readRDS("rma_HT_HG-U133A_EA.rds")
ht_hga <- readRDS("rma_HT_HG-U133A.rds")
hga <- readRDS("rma_HG-U133A.rds")

#log2 ht_hga (RMA from xps doesn't log2)
exprs(ht_hga) <- log2(exprs(ht_hga))

#fix up sample names
ht_hga_names <- strsplit(sampleNames(ht_hga), "[.]")
ht_hga_names <- sapply(ht_hga_names, function(x) substring(x[1], 2))
ht_hga_names <- gsub("_", ".", ht_hga_names)
sampleNames(ht_hga) <- ht_hga_names

sampleNames(ht_hga_ea) <- gsub(".CEL", "", sampleNames(ht_hga_ea))
sampleNames(hga) <- gsub(".CEL", "", sampleNames(hga))

#merge data from all platforms
all_exprs <- merge(exprs(hga), exprs(ht_hga), by="row.names")
row.names(all_exprs) <- all_exprs$Row.names
all_exprs <- merge(all_exprs, exprs(ht_hga_ea), by="row.names")
row.names(all_exprs) <- all_exprs$Row.names
all_exprs <- all_exprs[,-(1:2)]


#----------------
# SETUP ANALYSIS
#----------------

#generate model matrix
cmap_instances <- read.table("~/Documents/Batcave/GEO/2-cmap/data/raw/cmap_instances_02.csv",
                             header=T, sep="\t", quote='', fill=T, stringsAsFactors=F)

#"valid" drug names (needed for limma makeContrasts)
drugs <- unique(cmap_instances$cmap_name)
drugs_val <- make.names(drugs, unique=T)

samples <- colnames(all_exprs)
controls = c()

mod <- data.frame(matrix(data=0, nrow=length(samples), ncol=length(drugs)+1,
                         dimnames=list(samples, c(drugs_val, "ctl"))))

for (i in seq_along(drugs)) {
  drug <- drugs[i]
  drug_val <- drugs_val[i]

  #get cmap info for drug
  drug_instances <- cmap_instances[cmap_instances$cmap_name == drug, ]
  
  #add perturbation ids to model matrix
  pids <- drug_instances$perturbation_scan_id
  mod[pids, drug_val] <- 1
  
  #select control ids
  cids <- drug_instances$vehicle_scan_id4
  
  #generate full control id names
  for (i in seq_along(cids)) {
    cid <- cids[i]
    pid <- pids[i]
    
    #if multiple controls: get prefix/sufixes 
    if (length(strsplit(cid, "[.]")[[1]]) > 2) {
      
      pref = strsplit(pid,"[.]")[[1]][1]
      sufs = strsplit(cid,"[.]")
      sufs = sufs[[1]][-which(sufs[[1]] %in% "")]
      
      #paste prefix/suffixes together
      controls <- c(controls, paste(pref, sufs, sep="."))
      
      #if single control: add cid directly
    } else {
      controls <- unique(c(controls, cid))
    }
  }
}

#add ctls to mod
mod[controls, "ctl"] <-  1

#generate null model matrix for SVA
mod0 <- model.matrix(~1, data=mod)
row.names(mod) <- row.names(mod0)

#turn mod/all_exprs into matrix
all_exprs <- as.matrix(all_exprs)
mod <- as.matrix(mod)

#generate eset from all_exprs
eset <- new("ExpressionSet", exprs = all_exprs)


#-----------
# ANALYSIS
#-----------

#perform sva (2+ hours)
svobj <- sva(all_exprs, mod, mod0)

#add SVs to mod
modsv <- cbind(mod, svobj$sv)
colnames(modsv) <- c(colnames(mod), paste("SV", 1:svobj$n.sv, sep=""))

#generate contrast names (must be "valid")
contrasts <- paste(drugs_val, "ctl", sep="-")

#run limma analysis (2+ hours)
ebayes_sv <- fit_ebayes(eset, contrasts, modsv)


#---------
# COMBINE
#---------

#values to calc dprime
df <- ebayes_sv$df.residual + ebayes_sv$df.prior
ni <- sum(mod[, "ctl"])

top_tables <- data.frame(row.names=featureNames(eset))
for (i in seq_along(drugs)) {

  #get top table
  top_table <- topTable(ebayes_sv, coef=i, n=Inf)

  #add dprime
  nj <- sum(mod[, i])
  top_table$dprime <- effectsize(top_table$t, ((ni * nj)/(ni + nj)), df)[, "dprime"]

  #bind top tables (probe order from eset)
  colnames(top_table) <- paste(drugs[i], colnames(top_table), sep="_")
  top_tables <- cbind(top_tables, top_table[featureNames(eset),])
}

tops <- list(data=top_tables, drugs=drugs)
saveRDS(tops, "~/Documents/Batcave/GEO/2-cmap/data/processed/es/probes_top_tables.rds")


#collect dprime values (effect size) from each contrast
dp <- grepl("dprime|adj.P.Val", colnames(top_tables))
es_probes <- top_tables[, dp]
colnames(es_probes) <- gsub("_dprime", "", colnames(es_probes))

saveRDS(as.matrix(es_probes), 
        "~/Documents/Batcave/GEO/2-cmap/data/processed/es/probes_drug_es.rds")


#-----------
# ANNOTATE
#-----------

#get dprime and adj.P.Val from each contrast
dp <- grepl("dprime|adj.P.Val", colnames(top_tables))
es_probes <- top_tables[, dp]

#annotate results
library(hgu133a.db)
map <- AnnotationDbi::select(hgu133a.db, row.names(es_probes), "SYMBOL")
map <- map[!is.na(map$SYMBOL), ]
es_probes <- es_probes[map$PROBEID, ] #expands 1:many
es_probes[,"SYMBOL"] <- map$SYMBOL

#where SYMBOL duplicated, keep SYMBOL with smallest p-value
library(dplyr)
es_genes <- data.frame(SYMBOL=unique(map$SYMBOL), stringsAsFactors=F)

for (i in seq_along(drugs)) {

  #select adj.P.Val, dprime, and SYMBOL columns
  cols <- colnames(es_probes)[(i*2-1):(i*2)]
  es <- es_probes[, c(cols, "SYMBOL")]

  #wrap p-value column name in backticks (needed for arrange)
  ar <- paste("`", cols[1], "`", sep="")
  
  es %>%
    group_by(SYMBOL) %>%
    arrange_(.dots = ar) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    inner_join(es_genes, by="SYMBOL") ->
    es_genes
}

#symbol for row names
class(es_genes) <- "data.frame"
row.names(es_genes) <- es_genes$SYMBOL

#keep dprime columns only
dp <- grepl("dprime", colnames(es_genes))
es_genes <- es_genes[, dp]

#remove dprime from column names
colnames(es_genes) <- gsub("_dprime", "", colnames(es_genes))
es_genes <- es_genes[, drugs]
 
#save results
processed_rma <- list(eset=eset, svobj=svobj, ebayes_sv=ebayes_sv)

saveRDS(processed_rma, "processed_rma.rds")
saveRDS(as.matrix(es_genes), 
        "~/Documents/Batcave/GEO/2-cmap/data/processed/es/genes_drug_es.rds")