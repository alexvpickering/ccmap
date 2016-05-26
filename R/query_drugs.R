#' Get overlap between query and drug signatures.
#'
#' Determines the number of genes that change in the same direction between
#' query and drug signatures.
#'
#' Drugs with the largest positive and negative net overlap are predicted to,
#' respectively, mimic and reverse the query signature.
#'
#' @importFrom methods as
#' @importFrom utils head tail
#'
#' @param query_genes Named numeric vector of differentual expression values for
#'   query genes. Usually 'meta' slot of \code{get_dprimes} result.
#' @param drug_info Matrix of differential expression values or gene rankings
#'   (must set es=FALSE) for drugs. Gene rankings should be in increasing order
#'   (e.g. 1 is the most up or down regulated gene). Rows are genes, columns are
#'   drugs.
#' @param query_n Integer specifying how many of the top DE gene should be used
#'   from query_genes. If there are less than query_n genes in common between
#'   query and drug signatures, the number of common genes will be used.
#'   Parameter ignored if \code{step} not \code{NULL}.
#' @param drug_n Integer specifying how many of the top DE genes should be
#'   used from drug_genes. If there are less than query_n genes in common
#'   between query and drug signatures, the number of common genes will be used.
#' @param step Integer specifying step size for query range. Queries start
#'    at 100 genes, and step up by specified step size up to the number of
#'    common genes.
#' @param es Does drug_info contain effect size data? If FALSE, gene ranking
#'   data is assumed.
#' @family query functions
#'
#' @return List of data.frames (one per query size) each with columns:
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


query_drugs <- function(query_genes, drug_info = NULL,
                       query_n = length(query_genes),
                       drug_n = nrow(drug_info),
                       step = NULL,
                       es = TRUE) {

    #default to cmap_es for drug_info
    if (is.null(drug_info)) drug_info <- ccdata::cmap_es

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

#' Get overlap between query and drug signatures for a range of query gene sizes.
#'
#' Queries drug signatures using a range of query gene sizes. Results are sorted
#' by an auc metric which has the advantage of weighting query genes according
#' to their extent of differential expression.
#'
#' Drugs with, respectively, the largest positive and negative auc (head and
#' tail of resulting data.frame) are predicted to mimic and reverse the query
#' signature.
#'
#' @inheritParams query_drugs
#' @family query functions
#'
#' @return data.frame with number of net genes that overlap between query and
#'    drug signatures. Columns correspond to query sizes, rows to drugs. Drugs
#'    are sorted by decreasing area under the curve formed by plotting net
#'    overlap as a function of query size.
#'
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
#' res <- range_query_drugs(query_sig, as.matrix(drug_info))

range_query_drugs <- function(query_genes, drug_info=NULL, es=TRUE, step=100) {

    #query drugs for a range of query sizes
    top_drugs <- query_drugs(query_genes, drug_info, es=es, step=step)

    #create ranks and nets dataframes
    return(get_range_res(top_drugs))
}
