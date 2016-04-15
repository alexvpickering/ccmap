library(apcluster)
library(reshape2)
library(ggplot2)
library(data.table)
library(xgboost)

source("~/Documents/Batcave/GEO/2-cmap/combo_utils.R")

#---------------------

get_top_drugs <- function(query_genes, query_n=length(query_genes), drug_info=NULL,
                          drug_genes=NULL, query_type=NULL, drug_genes_n=NULL, es=F) {

  if (is.null(drug_genes)) {
    #get top up/dn drug genes
    drug_genes <- get_drug_genes(drug_info, query_genes, drug_genes_n, es)
  }
  #get up/dn query genes
  query_genes <- setup_query_genes(query_genes, query_n, drug_info, query_type)

  #determine direct overlap
  overlap <- get_overlap(drug_genes, query_genes)

  #get cross overlap (drug acts opposite to query)
  #SHOULD BE ABLE TO REMOVE AS CROSS IS DEFINED BY OVERLAP
  cross_genes <- list(up=query_genes$dn, dn=query_genes$up)
  cross_overlap <- get_overlap(drug_genes, cross_genes)

  #get net overlap
  overlap$table <- get_net_overlap(overlap, cross_overlap)

  return(overlap)
}

#---------------------

get_drug_genes <- function(drug_info, query_genes, drug_genes_n=NULL, es=F) {
  #used to get top genes_n up and down regulated genes in drug_info
  
  drug_genes <- list()
  drugs <- colnames(drug_info)

  #only consider drug genes that are also in query genes
  drug_info <- drug_info[row.names(drug_info) %in%
                         toupper(names(query_genes)), ]

  if (is.null(drug_genes_n)) {
    drug_genes_n <- nrow(drug_info)
  }

  for (drug in drugs) {


    if (es) {
      #get up/dn for each drug
      dr <- drug_info[, drug]
      updn <- get_updn(dr, drug_genes_n)

      up <- updn$up
      dn <- updn$dn
    } else {
      dr <- sort(drug_info[, drug])

      up <- names (head(dr, drug_genes_n / 2))
      dn <- names (tail(dr, drug_genes_n / 2))
    }
    
    #up/dn regulated
    
    top_list <- list(up=up, dn=dn)
    drug_genes[[drug]] <- top_list
  }
  return(drug_genes) 
}

#---------------------
get_updn <- function (query_genes, query_n) {
  #used by setup_query_genes to get up and dn

  #sort query genes
  query_genes <- sort(query_genes, decreasing=T)

  #get up/dn list
  up <- names(query_genes[query_genes > 0])
  dn <- names(query_genes[query_genes < 0])

  if (length(up) < query_n / 2 & 
      length(dn) > query_n / 2) {
      #take more genes from dn
      n_up <- length(up)
      n_dn <- query_n - length(up)
    } else if (length(dn) < query_n / 2 & 
               length(up) > query_n / 2) {
      #take more genes from up
      n_up <- query_n - length(dn)
      n_dn <- length(dn)
    } else {
      #take query_n / 2 from up and dn
      n_up <- query_n / 2
      n_dn <- query_n / 2
    }
  #get n_up, n_dn from up/dn list
  up <- head(up, n_up)
  dn <- tail(dn, n_dn)
  return(list(up=up, dn=dn))
}
#----------------

setup_query_genes <- function(query_genes, query_n=length(query_genes),
                              drug_info=NULL, query_type=NULL) {
  #used to format query genes appropriately

  #make names uppercase
  names(query_genes) <- toupper(names(query_genes))

  if (!is.null(drug_info)) {
    #remove genes not in drug_info
    query_genes <- query_genes[names(query_genes) %in% row.names(drug_info)]
  }

  if (is.null(query_type)) {

    updn <- get_updn(query_genes, query_n)
    up <- updn$up
    dn <- updn$dn  

  } else if (query_type == "prl") {

    n <- length(query_genes)
    up <- names (head(query_genes[1:n/2], query_n / 2))
    dn <- names (tail(query_genes[n/2:n], query_n / 2))
  }
  avail_n <- length(up) + length(dn)

  query_list <- list(up=up, dn=dn)
  return(query_list)
}

#---------------------

get_overlap <- function(drug_genes, query_genes, tan_sim=FALSE) {
  #used to provide overlap between drug_genes and query_genes

  drugs <- names(drug_genes)
  overlap_table <- data.frame(row.names=drugs, stringsAsFactors=F)
  overlap_genes <- list()

  for (dr in drugs) {
    dr_genes <- drug_genes[[dr]]

    #get up/down and tot overlap
    up <- intersect(dr_genes$up, query_genes$up)
    dn <- intersect(dr_genes$dn, query_genes$dn)
    tot <- length(up) + length(dn)

    #squirel away
    overlap <- list(up=up, dn=dn)
    overlap_genes[[dr]] <- overlap

    if (tan_sim) {
      #overlap is tan coeficient
      unique <- length (union(unlist(dr_genes), unlist(query_genes)))
      overlap_table[dr, "overlap"] <- tot / unique
    } else {
      #overlap is total count
      overlap_table[dr, "overlap"] <- tot

      #add number of drug/query genes
      drug_genes_n <- length(dr_genes$up) + length(dr_genes$dn)
      query_genes_n <- length(query_genes$up) + length(query_genes$dn)

      overlap_table[dr, "drug_genes_n"] <- drug_genes_n
      overlap_table[dr, "query_genes_n"] <- query_genes_n
    }
  }
  overlap <- list(table=overlap_table, genes=overlap_genes)
  return (overlap)
}

#---------------------

get_net_overlap <- function(overlap, cross, tan_sim=F) {
  #returns overlap with table sorted by net overlap
  table <- overlap$table

  #add cross overlap to table
  table[, "cross"] <- cross$table$overlap

  #compute net overlap: overlap - cross
  net <- table$overlap - table$cross
  table[, "net"] <- net

  if (!tan_sim) {
    #add net_query_pct and sort by net
    pct <- round(net / table$query_genes_n * 100, 2)
    table[, "net_query_pct"] <- pct
    table <- table[order(net, decreasing=T), ]
    #reorder columns
    table <- table[,c(2,1,4,5,3,6)]
  } 
  return(table)
}

#---------------------


get_dprimes <- function(es2) {
  #used to extract meta/contrast dprimes from ES.GeneMeta result
    
  scores_table <- es2$theScores

  meta <- scores_table[, "MUvals"]
  contrasts <- list()

  ex_cols <- grep("Effect_Ex_", colnames(scores_table))
  
  for (col in ex_cols){
      
      dprimes <- scores_table[, col]
      col_name <- colnames(scores_table)[col]
      contrasts[[col_name]] <- dprimes
  }
  return(list (meta=meta, contrasts=contrasts))  
}

#---------------------

get_tstats <- function(es, diff_exprs) {

  #get meta t-stats
  meta <- es$TestStatistic
  names(meta) <- featureNames(diff_exprs[[1]]$eset)

  #get contrast t-stats
  contrasts <- list()
  top_tables <- get_tts(diff_exprs)

  for (i in seq_along(top_tables)) {

    tstats <- top_tables[[i]]$t
    names(tstats) <- row.names(top_tables[[i]])
    contrasts[[i]] <- tstats
  }
  return(list (meta=meta, contrasts=contrasts))  
}

#---------------------

get_prlscores <- function(prl, diff_exprs) {
  #uses tstats for contrast scores (used by merge_ranks to obtain prl)
  
  meta <- prl
  
  #get contrast t-stats
  contrasts <- list()
  top_tables <- get_tts(diff_exprs)

  for (i in seq_along(top_tables)) {

    tstats <- top_tables[[i]]$t
    names(tstats) <- row.names(top_tables[[i]])
    contrasts[[i]] <- tstats
  }
  return(list (meta=meta, contrasts=contrasts))  
}

#---------------------

test_ma <- function(scores, drug_name, drug_info, query_ns=2^(1:13),
                    query_type=NULL, drug_genes_n=NULL, es=F) {

  #get top up/dn drug genes
  drug_genes <- get_drug_genes(drug_info, query_genes, drug_genes_n, es)

  #final results 
  res_cons <- list()
  res_meta <- c()

  for (query_n in query_ns) {

    #perform query of size query_n with meta result
    meta_top <- get_top_drugs(scores$meta, query_n, drug_info,
                              drug_genes, query_type, drug_genes_n, es)

    meta_rank <- get_drug_rank(meta_top, drug_name)
    res_meta <- c(res_meta, meta_rank)

    cons_ranks <- c()
    for (con in scores$contrasts) {

      #perform query of size query_n with each contrast
      con_top <- get_top_drugs(con, query_n, drug_info, 
                               drug_genes, query_type, drug_genes_n, es)

      con_rank <- get_drug_rank(con_top, drug_name)

      #combine rank with those of other contrasts
      cons_ranks <- c(cons_ranks, con_rank)
    }
    #store ranks of each contrast for size query_n
    res_cons[[as.character(query_n)]] <- cons_ranks
  }
  #add names to res_meta
  names(res_meta) <- as.character(query_ns)

  ma_res <- list(meta=res_meta, cons=res_cons)
  return(ma_res)
}

#---------------------

get_drug_rank <- function(top_drugs, drug_name) {
  #used by test_ma: returns number of drugs with net >= drug_name

  table <- top_drugs$table

  drug_net <- table[drug_name, "net"]
  drug_rank <- nrow(table[table$net >= drug_net, ])

  return(drug_rank)
}

#---------------------

plot_ma_res <- function (ma_res) {

  cons_res <- data.frame(ma_res$cons)
  meta_res <- t(data.frame(ma_res$meta))
  colnames(cons_res) <- colnames(meta_res)

  cons_res <- melt(cons_res)
  meta_res <- melt(meta_res)[,-1]
  colnames(cons_res) <- c("query_n", "rank")
  colnames(meta_res) <- c("query_n", "rank")

  meta_res$meta <- 1
  cons_res$meta <- 0

  full_res <- rbind(meta_res, cons_res)
  full_res$rank <- as.integer(full_res$rank)
  full_res$query_n <- as.integer(full_res$query_n)

  p <- ggplot(full_res, aes(log2(query_n), rank))
  p +
    stat_summary(fun.y="median", colour = "red", size = 10, geom="point", shape=45) +
    geom_point(aes(alpha=factor(meta))) + 
    scale_y_continuous(breaks = seq(0, 1310, 100), limits=c(0, 1310)) +
    scale_x_continuous(breaks = seq(1, 13, 1), limits=c(1, 13))

}


#---------------------

get_top_combos <- function(query_genes, query_n=length(query_genes),
                           drug_info=NULL, drug_genes_n=NULL) {

  #drug_info: ONLY drug_es (not drug_prls) 

  #--------------------
  # ALGORITHM (greedy):
  #--------------------

  #  - get best drug
  #     > combine with each other drug (probe level)
  #     > convert probes to genes
  #  - get best combo

  #--------
  # SETUP
  #--------

  #load model & cmap data
  mod    <- readRDS("~/Documents/Batcave/GEO/2-cmap/data/combos/model.rds")
  cmap   <- readRDS("~/Documents/Batcave/GEO/2-cmap/data/processed/es/probes_top_tables.rds")
  data   <- cmap$data
  drugs  <- cmap$drugs
  probes <- row.names(data)

  #indices for drugs
  inds <- sapply(seq_along(drugs), function(x) c(x*7-6, x*7))
  colnames(inds) <- drugs

  #start combo table with top drug
  combo_table <- get_top_drugs(query_genes, query_n, drug_info,
                               drug_genes_n=drug_genes_n, es=T)$table[1, ]

  top_drug <- row.names(combo_table)
  other_drugs <- setdiff(drugs, top_drug)

  #get combo data
  pairs <- lapply(seq_along(other_drugs), function(x) c(top_drug, other_drugs[x]))
  names(pairs) <- sapply(pairs, function(x) paste(x[1], x[2], sep=" + "))
  combo_data <- combine_cmap_pairs(pairs, data, inds)

  #load annotation information
  suppressMessages(library("hgu133a.db"))
  suppressMessages(map <- AnnotationDbi::select(hgu133a.db, probes, "SYMBOL"))
  map <- map[!is.na(map$SYMBOL),]
  genes <- unique(map$SYMBOL)

  #--------------
  # PROBE PREDS
  #--------------

  #make predictions with mod
  probe_preds <- predict(mod, combo_data) - 0.5

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

  #where duplicated SYMBOL, choose probe with prediction furthest from 0.5
  probe_preds <- data.table(probe_preds)
  gene_preds  <- probe_preds[, lapply(.SD, function (col) col[which.max(abs(col))]), by='SYMBOL']

  #format table
  gene_preds <- as.data.frame(gene_preds)
  row.names(gene_preds) <- gene_preds$SYMBOL
  gene_preds <- gene_preds[, colnames(gene_preds) != "SYMBOL"]

  #------------
  # BEST COMBO
  #------------

  #get best combo of top_drugs and other_drugs
  best_combos <- get_top_drugs(query_genes, query_n, drug_info=as.matrix(gene_preds),
                               drug_genes_n=drug_genes_n, es=T)$table

  combo_table <- rbind(combo_table, best_combos)
  return(combo_table)
}