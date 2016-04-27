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



get_range_res <- function(test_res) {

    #setup ranks and nets dataframes
    drugs <- row.names(test_res[[1]])
    query_ns <- as.integer(names(test_res))

    ranks <- data.frame(row.names=drugs)
    nets  <- data.frame(row.names=drugs)

    #fill in dataframes
    for (i in seq_along(test_res)) {

        query_n <- names(test_res)[i]
        ranks[drugs, query_n] <- sapply(drugs, get_drug_rank, test_res[[i]])
        nets[drugs, query_n] <- test_res[[i]][drugs, "net"]
    }


    #sort nets by auc formed by plot of nets vs query_n
    aucs <- sapply(as.data.frame(t(nets)), function(y) MESS::auc(query_ns, y))
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





clean_rbindlist <- function(l) {

    #get row names
    l_rnames <- unlist(lapply(l, row.names))

    #bind dfs
    df <- data.table::rbindlist(l)

    #clean up
    class(df) <- "data.frame"
    row.names(df) <- l_rnames

    return(df)
}

