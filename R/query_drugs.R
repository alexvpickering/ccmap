#' Get overlap between query and drug signatures.
#'
#' Determines the number of genes that change in the same direction between
#' query and drug signatures.
#'
#' Drugs with the largest positive and negative net overlap are predicted to,
#' respectively, mimic and reverse the query signature.
#'
#' @param query_genes Named numeric vector of differentual expression values for
#'   query genes. Usually 'meta' slot of \code{get_dprimes} result.
#' @param drug_info Matrix of differential expression values for drugs or drug
#'   combinations. Rows are genes, columns are drugs.
#' @param query_n Integer specifying how many of the top DE gene should be used
#'   from query_genes. If there are less than query_n genes in common between
#'   query and drug signatures, the number of common genes will be used.
#'   Parameter ignored if step not NULL.
#' @param drug_n Integer specifying how many of the top DE genes should be
#'   used from drug_genes. If there are less than query_n genes in common
#'   between query and drug signatures, the number of common genes will be used.
#' @param step Integer specifying step size for query range. Queries start
#'    at 100 genes, and step up by specified step size up to the number of
#'    common genes.
#'
#' @seealso \code{\link{predict_combos}} to get predicted drug combination
#'   signatures (can be passed to \code{drug_info}).
#'
#' @return If step is not NULL, a data.frame sorted by auc of net overlap vs
#'    query size. Columns correspond to query sizes, rows to drugs.
#'
#'    If step is NULL, a data.frame sorted by net overlap with columns:
#'   \item{drug_n}{Number of drug genes used.}
#'   \item{query_n}{Number of query genes used.}
#'   \item{overlap}{Number of genes that change in the same direction in drug
#'      and query signatures.}
#'   \item{cross}{Number of genes that change in the opposite direction in drug
#'      and query signatures.}
#'   \item{net}{Difference between overlap and cross.}
#' @export
#'
#' @examples
#'
#' # create drug signatures
#' genes <- paste("GENE", 1:1000, sep = "_")
#' set.seed(0)
#'
#' drug_info <- data.frame(row.names = genes,
#'                         drug1  = rnorm(1000, sd = 2),
#'                         drug2  = rnorm(1000, sd = 2),
#'                         drug3  = rnorm(1000, sd = 2))
#'
#' # query signature is drug3
#' query_sig <- drug_info$drug3
#' names(query_sig) <- genes
#'
#' res <- query_drugs(query_sig, as.matrix(drug_info))


query_drugs <- function(query_genes, drug_info = NULL, step = NULL,
                        query_n = length(query_genes), drug_n = nrow(drug_info)) {
    #bind global
    cmap_es = NULL

    #default to cmap_es for drug_info
    if (is.null(drug_info)) {
        utils::data("cmap_es", package = "ccdata", envir = environment())
        drug_info <- cmap_es
    }

    #get top up/dn drug genes
    drug_genes <- get_drug_genes(drug_info, query_genes, drug_n)

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

    if (is.null(step)) return(top_drugs[[1]])

    return(get_range_res(top_drugs))
}
