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

query_drugs <- function(query_genes, drug_info,
                       query_n=length(query_genes),
                       drug_n=nrow(drug_info),
                       step=NULL,
                       es=TRUE) {

    #get top up/dn drug genes
    drug_genes <- get_drug_genes(drug_info, query_genes, drug_n, es)

    if (!is.null(step)) {
        #setup range of query sizes
        max_n <- length(unlist(drug_genes[[1]]))
        query_n <- round(seq(100, round(max_n, -2), step))
        query_n[length(query_n)] <- max_n
    }

    #get overlap df for each query_n
    top_drugs <- list()

    for (n in query_n) {
        #get up/dn query genes
        q_genes <- setup_query_genes(query_genes, n, drug_info)

        #determine direct overlap
        overlap <- get_overlap(drug_genes, q_genes)

        #calculate cross and net overlap
        overlap$cross <- overlap$query_n - overlap$overlap
        overlap$net   <- overlap$overlap - overlap$cross

        #sort rows and reorder columns
        overlap <- overlap[order(overlap$net, decreasing=TRUE),
                           c("drug_n", "query_n", "overlap", "cross", "net")]

        used_n <- as.character(overlap$query_n[1])
        top_drugs[[used_n]] <- overlap
    }
    return(top_drugs)
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

range_query_drugs <- function(query_genes, drug_info, es=TRUE, step=100) {

    #query drugs for a range of query sizes
    top_drugs <- query_drugs(query_genes, drug_info, es=es, step=step)

    #create ranks and nets dataframes
    return(get_range_res(top_drugs))
}
