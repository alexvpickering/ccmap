#' Get overlap between query and drug signatures.
#'
#'
#'
#' @param query_genes Named numeric vector of differentual expression values for
#'   query genes. Usually 'meta' slot of \code{get_dprimes} result.
#' @param query_n Integer specifying how many of the top DE gene should be used
#'   from query_genes. If there are less than query_n genes in common between
#'   query and drug signatures, the number of common genes will be used.
#' @param drug_info Matrix of differential expression values or gene rankings
#'   (must set es=FALSE) for drugs. Gene rankings should be in increasing order
#'   (e.g. 1 is the most up or down regulated gene). Rows are genes, columns are
#'   drugs.
#' @param drug_genes List of lists (one per drug) with slots 'up' and 'dn' each
#'   containing a character vector with the names of up- and down-regulated
#'   genes. If not supplied, created by \code{get_drug_genes}. Useful only if
#'   multiple calls to \code{get_top_drugs} needed (e.g. in \code{test_ma}).
#' @param drug_n Integer specifying how many of the top DE genes should be
#'   used from drug_genes. If there are less than query_n genes in common
#'   between query and drug signatures, the number of common genes will be used.
#' @param es Does drug_info contain effect size data? If FALSE, gene ranking
#'   data is assumed.
#'
#' @return data.frame with columns:
#'   \item{}{}
#'   \item{}{}
#'   \item{}{}
#'   \item{}{}
#'   \item{}{}
#' @export
#'
#' @examples

get_top_drugs <- function(query_genes, drug_info,
                          query_n=length(query_genes), drug_genes=NULL,
                          drug_n=nrow(drug_info), es=TRUE) {

    if (is.null(drug_genes)) {
        #get top up/dn drug genes
        drug_genes <- get_drug_genes(drug_info, query_genes, drug_n, es)
    }
    #get up/dn query genes
    query_genes <- setup_query_genes(query_genes, query_n, drug_info)

    #determine direct overlap
    overlap <- get_overlap(drug_genes, query_genes)

    #calculate cross and net overlap
    overlap$cross <- overlap$query_n - overlap$overlap
    overlap$net   <- overlap$overlap - overlap$cross

    #sort rows and reorder columns
    overlap <- overlap[order(overlap$net, decreasing=T),
                       c("drug_n", "query_n", "overlap", "cross", "net")]

    return(overlap)
}

#---------------------

#' Get list of top differentially regulated genes for drugs.
#'
#' @inheritParams get_top_drugs
#'
#' @export
#' @seealso get_top_drugs
#' @return List of lists (one per drug in drug_info) with slots:
#'   \item{up}{Named vector of genes up-regulated by drug.}
#'   \item{up}{Named vector of genes down-regulated by drug.}
#'

get_drug_genes <- function(drug_info, query_genes,
                           drug_n=nrow(drug_info), es=TRUE) {

    drug_genes <- list()
    drugs <- colnames(drug_info)

    #only consider drug genes that are also in query genes
    drug_info <- drug_info[row.names(drug_info) %in%
                               toupper(names(query_genes)), ]

    #max drug_n is drug genes that are also in query genes
    if (drug_n > nrow(drug_info)) {
        drug_n <- nrow(drug_info)
    }
    for (drug in drugs) {
        if (es) {
            #get up/dn for each drug
            dr <- drug_info[, drug]
            updn <- get_updn(dr, drug_n)

            up <- updn$up
            dn <- updn$dn
        } else {
            dr <- sort(drug_info[, drug])

            up <- names (head(dr, drug_n / 2))
            dn <- names (tail(dr, drug_n / 2))
        }
        #up/dn regulated
        top_list <- list(up=up, dn=dn)
        drug_genes[[drug]] <- top_list
    }
    return(drug_genes)
}

#----------------

#' Get top DE genes for query.
#'
#' @inheritParams get_top_drugs
#'
#' @return List with slots 'up' and 'dn' each containing captalized character
#'   vectors of up- and down-regulated genes.

setup_query_genes <- function(query_genes, query_n, drug_info) {

    #make names uppercase
    names(query_genes) <- toupper(names(query_genes))

    #remove genes not in drug_info
    query_genes <- query_genes[names(query_genes) %in% row.names(drug_info)]

    updn <- get_updn(query_genes, query_n)

    query_list <- list(up=updn$up, dn=updn$dn)
    return(query_list)
}

#---------------------

#' Get number of top DE genes.
#'
#' Function returns at most the requested number of genes, attempting to devide
#' equally between the top up and down regulated genes.
#'
#' If fewer than half of requested genes are available among up (or down)
#' regulated genes, then more genes will be taken from down (or up) regulated
#' genes.
#'
#' @param genes Named numeric vector of gene effect sizes.
#' @param genes_n Integer specifying the desired number of the top DE genes.
#'
#' @seealso \code{\link{setup_query_genes}}, \code{\link{get_drug_genes}}
#' @return List with slots 'up' and 'dn' each containing character vectors
#'   of up- and down-regulated genes.

get_updn <- function (genes, genes_n) {

    #sort query genes
    genes <- sort(genes, decreasing=TRUE)

    #get up/dn list
    up <- names(genes[genes > 0])
    dn <- names(genes[genes < 0])

    if (length(up) < genes_n / 2 &
        length(dn) > genes_n / 2) {
        #take more genes from dn
        n_up <- length(up)
        n_dn <- genes_n - length(up)
    } else if (length(dn) < genes_n / 2 &
               length(up) > genes_n / 2) {
        #take more genes from up
        n_up <- genes_n - length(dn)
        n_dn <- length(dn)
    } else {
        #take genes_n / 2 from up and dn
        n_up <- genes_n / 2
        n_dn <- genes_n / 2
    }
    #get n_up, n_dn from up/dn list
    up <- head(up, n_up)
    dn <- tail(dn, n_dn)

    return(list(up=up, dn=dn))
}

#---------------------

#' Title
#'
#' @param drug_genes
#' @param query_genes
#' @param tan_sim
#'
#' @return
#' @export
#'
#' @examples
get_overlap <- function(drug_genes, query_genes, tan_sim=FALSE) {
    #used to provide overlap between drug_genes and query_genes

    drugs <- names(drug_genes)
    overlap <- data.frame(row.names=drugs, stringsAsFactors=FALSE)

    #add number of drug/query genes
    overlap$query_n <- length(query_genes$up) + length(query_genes$dn)
    overlap$drug_n  <- length(drug_genes[[1]]$up) + length(drug_genes[[1]]$dn)

    for (dr in drugs) {
        dr_genes <- drug_genes[[dr]]

        #get up/down and total overlap
        up <- intersect(dr_genes$up, query_genes$up)
        dn <- intersect(dr_genes$dn, query_genes$dn)
        tot <- length(up) + length(dn)

        #add total overlap for drug
        overlap[dr, "overlap"] <- tot

    }
    return (overlap)
}

#---------------------


#' Extract unbiased effect sizes from meta-analysis by ES.GeneMeta.
#'
#' Function extracts MUvals (overall mean effect size) and Effect_Ex (unbiased
#' effect sizes from each contrast).
#'
#' Result used to query connectivity map drugs.
#'
#' @param es2 ES.GeneMeta.res object, result of of call to \code{ES.GeneMeta}
#'
#' @return List containing:
#'   \item{meta}{Named numeric vector with overall mean effect size for each gene
#'      from meta analysis.}
#'   \item{contrasts}{List of named numeric vectors (one per contrast) with
#'      unbiased effect sizes for each gene from meta analysis.}
#' @export
#'
#' @seealso \link{\code{ES.GeneMeta}}.
#' @examples
get_dprimes <- function(es2) {

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


#' Title
#'
#' @param scores
#' @param drug_name
#' @param drug_info
#' @param query_ns
#' @param es
#'
#' @return
#' @export
#'
#' @examples
test_ma <- function(scores, drug_name, drug_info, query_ns=2^(4:13), es=TRUE) {

    #get top up/dn drug genes
    drug_genes <- get_drug_genes(drug_info, scores$meta, es=es)

    #final results
    res_cons <- list()
    res_meta <- c()

    for (query_n in query_ns) {

        #perform query of size query_n with meta result
        meta_top <- get_top_drugs(scores$meta, drug_info,
                                  query_n, drug_genes, es=es)

        meta_rank <- get_drug_rank(drug_name, meta_top)
        res_meta <- c(res_meta, meta_rank)

        cons_ranks <- c()
        for (con in scores$contrasts) {

            #perform query of size query_n with each contrast
            con_top <- get_top_drugs(con, drug_info,
                                     query_n, drug_genes, es=es)

            con_rank <- get_drug_rank(drug_name, con_top)

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

#' Title
#'
#' @param top_drugs
#' @param drug_name
#'
#' @return
#' @export
#'
#' @examples
get_drug_rank <- function(drug_name, top_drugs) {
    #used by test_ma: returns number of drugs with net >= drug_name

    drug_net <- top_drugs[drug_name, "net"]
    drug_rank <- nrow(top_drugs[top_drugs$net >= drug_net, ])

    return(drug_rank)
}

#---------------------

#' Title
#'
#' @import ggplot2
#'
#' @param ma_res
#'
#' @return
#' @export
#'
#' @examples
plot_ma_res <- function (ma_res, ymax=1310, ytick=100) {

    cons_res <- data.frame(ma_res$cons)
    meta_res <- t(data.frame(ma_res$meta))
    colnames(cons_res) <- colnames(meta_res)

    cons_res <- reshape2::melt(cons_res)
    meta_res <- reshape2::melt(meta_res)[,-1]
    colnames(cons_res) <- c("query_n", "rank")
    colnames(meta_res) <- c("query_n", "rank")

    meta_res$meta <- 1
    cons_res$meta <- 0

    full_res <- rbind(meta_res, cons_res)
    full_res$rank <- as.integer(full_res$rank)
    full_res$query_n <- as.integer(full_res$query_n)

    p <- ggplot(full_res, aes(query_n, rank))

    p + stat_summary(fun.y="median", colour = "red",
                     size=10, geom="point", shape=45) +

        geom_jitter(width=0.3, aes(alpha=factor(meta))) +

        scale_y_continuous(breaks=seq(0, 1310, ytick),
                           limits=c(0, 1310)) +

        scale_x_continuous(breaks=2^(4:13),
                           limits=c(12, 10000),
                           trans=scales::log2_trans()) +

        coord_cartesian(ylim=c(0, ymax))
}

#---------------------

#' Title
#'
#' @param scores
#' @param drug_info
#' @param query_ns
#' @param es
#'
#' @return
#' @export
#'
#' @examples
test_ma_range <- function(dprimes, drug_info, es=TRUE) {

    drugs <- colnames(drug_info)

    #get top up/dn drug genes
    drug_genes <- get_drug_genes(drug_info, dprimes$meta, es=es)

    #query using genes numbers from 100 to max
    mx <- length(unlist(drug_genes[[1]]))
    query_ns <- round(seq(100, round(mx, -2), 100))
    query_ns[length(query_ns)] <- mx

    ranks <- data.frame(row.names=drugs)
    nets  <- data.frame(row.names=drugs)

    for (query_n in query_ns) {
        #perform query of size query_n with meta result
        top_drugs <- get_top_drugs(dprimes$meta, drug_info,
                                   query_n, drug_genes, es=es)

        #add meta rank and net for each drug
        col <- as.character(query_n)
        ranks[drugs, col] <- sapply(drugs, get_drug_rank, top_drugs)
        nets[drugs, col] <- top_drugs[drugs, "net"]
    }

    #sort nets by auc formed by plot of nets vs query_n
    aucs <- sapply(as.data.frame(t(nets)),
                   function(y) MESS::auc(query_ns, y))

    aucs <- sort(aucs, decreasing=TRUE)
    nets <- nets[names(aucs), ]


    #sort ranks by frequency of top ten ranks
    top_tens <- rowSums(ranks <= 10)
    top_tens <- sort(top_tens, decreasing=TRUE)

    ranks <- ranks[names(top_tens), ]

    return(list(ranks=ranks, nets=nets))
}




#---------------------

#' Title
#'
#' @import data.table
#' @import hgu133a.db
#'
#' @param query_genes
#' @param query_n
#' @param drug_info
#' @param drug_n
#'
#' @return
#' @export
#'
#' @examples
get_top_combos <- function(query_genes, query_n=length(query_genes),
                           drug_info=NULL, drug_n=NULL) {

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
                                 drug_n=drug_n)[1, ]

    top_drug <- row.names(combo_table)
    other_drugs <- setdiff(drugs, top_drug)

    #get combo data
    pairs <- lapply(seq_along(other_drugs), function(x) c(top_drug, other_drugs[x]))
    names(pairs) <- sapply(pairs, function(x) paste(x[1], x[2], sep=" + "))
    combo_data <- combine_cmap_pairs(pairs, data, inds)

    #load annotation information
    suppressMessages(map <- AnnotationDbi::select(hgu133a.db, probes, "SYMBOL"))
    map <- map[!is.na(map$SYMBOL),]
    genes <- unique(map$SYMBOL)

    #--------------
    # PROBE PREDS
    #--------------

    #make predictions with mod
    probe_preds <- xgboost::predict(mod, combo_data) - 0.5

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
    gene_preds <- probe_preds[,
                              lapply(.SD, function (col) col[which.max(abs(col))]),
                              by='SYMBOL']

    #format table
    gene_preds <- as.data.frame(gene_preds)
    row.names(gene_preds) <- gene_preds$SYMBOL
    gene_preds <- gene_preds[, colnames(gene_preds) != "SYMBOL"]

    #------------
    # BEST COMBO
    #------------

    #get best combo of top_drugs and other_drugs
    best_combos <- get_top_drugs(query_genes, query_n, drug_info=as.matrix(gene_preds),
                                 drug_n=drug_n)

    combo_table <- rbind(combo_table, best_combos)
    return(combo_table)
}
