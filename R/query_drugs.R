#' Get overlap between query and drug signatures.
#'
#' Determines the cosine similarity between the query and each drug signature.
#'
#' Drugs with the largest positive and negative cosine similarity are predicted to,
#' respectively, mimic and reverse the query signature. Values range from +1 to -1.
#'
#' The 230829 LINCS l1000 signatures (drugs & genetic over/under expression) can also be queried.
#' In order to compare l1000 results to those obtained with cmap, only the same genes should be
#' included (see example).
#'
#' @import ccdata
#' @param query_genes Named numeric vector of differentual expression values for
#'   query genes. Usually 'meta' slot of \code{get_dprimes} result.
#' @param drug_info Character vector specifying which dataset to query
#'   (either 'cmap' or 'l1000'). Can also provide a matrix of differential expression
#'   values for drugs or drug combinations (rows are genes, columns are drugs).
#' @param sorted Would you like the results sorted by decreasing similarity?
#'   Default is TRUE.
#'
#' @seealso \code{\link{query_combos}} to get similarity between query and
#'   predicted drug combination signatures.
#'
#' @return Vector of cosine similarities between query and drug combination signatures.
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
#' res <- query_drugs(query_sig, as.matrix(drug_info))
#'
#' # use only common genes for l1000 and cmap matrices
#' # library(ccdata)
#' # data(cmap_es)
#' # data(l1000_es)
#' # cmap_es <- cmap_es[row.names(l1000_es), ]


query_drugs <- function(query_genes, drug_info = c('cmap', 'l1000'), sorted = TRUE, ngenes = 100) {

    # default to cmap_es for drug_info
    if (class(drug_info) == 'character') {
        fname <- paste0(drug_info[1], '_es')
        utils::data(list = fname, package = "ccdata", envir = environment())
        drug_info <- get(fname)
        rm(list = fname)
    }

    # use only common genes
    query_genes <- query_genes[names(query_genes) %in% row.names(drug_info)]

    # top 100 up/down genes
    query_genes <- sort(query_genes, TRUE)
    query_genes <- c(head(query_genes, ngenes), tail(query_genes, ngenes))
    drug_info   <- drug_info[names(query_genes), ,drop = FALSE]

    # cosine similarity
    sim <- apply(drug_info, 2, lsa::cosine, query_genes)

    if (sorted) {
        return(sort(sim, decreasing = TRUE))
    } else {
        return(sim)
    }

}


#--------------------------------


#' Sum of cumulative sum computed over rows then columns of matrix.
#'
#' Equivalent to computing the cumulative sum of a matrix over rows, then
#' over columns, then suming every value (though much faster and more memory
#' efficient).
#'
#' @param x Numeric vector of non-zero values of matrix.
#' @param i Integer vector of row indices of x.
#' @param j Integer vector of column indices of x.
#'
#' @return Numeric value equal to the sum of the cumulative sum computed over
#'    rows then columns of a matrix.
#' @export
#'
#' @examples
#' x <- c(1, 1, 1, -1) # non-zero values of matrix
#' i <- c(1, 2, 3, 4)  # row indices of x
#' j <- c(4, 1, 3, 2)  # col indices of x
#'
#' sum_rowcolCumsum(x, i, j)

sum_rowcolCumsum <- function(x, i, j) {
    sum(x * (max(i) - i + 1) * (max(j) - j + 1))
}
