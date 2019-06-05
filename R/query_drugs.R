#' Get correlation between query and drug signatures.
#'
#' Determines the pearson correlation between the query and each drug signature.
#'
#' Drugs with the largest positive and negative pearson correlation are predicted to,
#' respectively, mimic and reverse the query signature. Values range from +1 to -1.
#'
#' The 230829 LINCS l1000 signatures (drugs & genetic over/under expression) can also be queried.
#' In order to compare l1000 results to those obtained with cmap, only the same genes should be
#' included (see second example).
#'
#' @import ccdata
#' @param query_genes Named numeric vector of differentual expression values for
#'   query genes. Usually 'meta' slot of \code{get_dprimes} result.
#' @param drug_info Character vector specifying which dataset to query
#'   (either 'cmap' or 'l1000'). Can also provide a matrix of differential expression
#'   values for drugs or drug combinations (rows are genes, columns are drugs).
#' @param sorted Would you like the results sorted by decreasing similarity?
#'   Default is TRUE.
#' @param ngenes The number of top differentially-regulated (up and down) query genes
#'   to use if \code{path} is NULL. If \code{path} is not NULL, \code{ngenes} is the larger
#'   of 15 or the number of pathway genes with absolute dprimes > 2.
#' @param path Character vector specifying KEGG pathway. Used to find drugs that most
#'   closely mimic or reverse query signature for specific pathway.
#'
#' @seealso \code{\link{query_combos}} to get similarity between query and
#'   predicted drug combination signatures. \link[crossmeta]{diff_path} and \link[crossmeta]{path_meta}
#'   to perform pathway meta-analysis.
#'
#' @return Vector of pearson correlations between query and drug combination signatures.
#'
#' @export
#'
#' @examples
#'
#' # Example 1 -----
#'
#' library(crossmeta)
#' library(ccdata)
#' library(lydata)
#'
#' data_dir <- system.file("extdata", package = "lydata")
#' data(cmap_es)
#'
#' # gather GSE names
#' gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")
#'
#' # load previous differential expression analysis
#' anals <- load_diff(gse_names, data_dir)
#'
#' # run meta-analysis
#' es <- es_meta(anals)
#'
#' # get meta-analysis effect size values
#' dprimes <- get_dprimes(es)
#'
#' # most significant pathway (from path_meta)
#' path <- 'Amino sugar and nucleotide sugar metabolism'
#'
#' # query using entire transcriptional profile
#' topd <- query_drugs(dprimes$all$meta, cmap_es)
#'
#' # query restricted to transcriptional profile for above pathway
#' topd_path <- query_drugs(dprimes$all$meta, cmap_es, path=path)
#'
#' # Example 2 -----
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


query_drugs <- function(query_genes, drug_info = c('cmap', 'l1000'), sorted = TRUE, ngenes = 200, path = NULL) {


    # bindings to pass check
    gslist = gs.names = NULL

    # default to cmap_es for drug_info
    if (class(drug_info) == 'character') {
        fname <- paste0(drug_info[1], '_es')
        utils::data(list = fname, package = "ccdata", envir = environment())
        drug_info <- get(fname)
        rm(list = fname)
    }

    if (!is.null(path)) {
        utils::data("gslist", "gs.names", package = "crossmeta", envir = environment())

        # path symbols
        path_num <- names(gs.names)[gs.names == path]
        path_sym <- unique(names(gslist[[path_num]]))

        # drug info with only path symbols
        drug_info <- drug_info[row.names(drug_info) %in% path_sym, ]
    }

    # use only common genes
    query_genes <- query_genes[names(query_genes) %in% row.names(drug_info)]

    # top up/down ngenes
    top_ngenes  <- utils::head(names(sort(abs(query_genes), TRUE)), ngenes)
    query_genes <- query_genes[top_ngenes]
    drug_info   <- drug_info[names(query_genes), ,drop = FALSE]

    # pearson correlation
    sim <- stats::cor(query_genes, drug_info, method="pearson")
    sim <- structure(c(sim), names=colnames(sim))

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
