#' Extract unbiased effect sizes from meta-analysis by crossmeta.
#'
#' Function extracts mu (overall mean effect size) and dprimes (unbiased
#' effect sizes from each contrast).
#'
#' Result used to query connectivity map drugs and predicted drug combinations.
#'
#' @param es Result of call to \code{es_meta}.
#'
#' @return List containing:
#'   \item{meta}{Named numeric vector with overall mean effect sizes for all genes
#'      from meta-analysis.}
#'   \item{contrasts}{List of named numeric vectors (one per contrast) with
#'      unbiased effect sizes for all measured genes.}
#' @export
#'
#' @seealso \code{\link[crossmeta]{es_meta}}.
#' @examples
#' library(crossmeta)
#' library(lydata)
#'
#' data_dir <- system.file("extdata", package = "lydata")
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
#' #get dprimes
#' dprimes <- get_dprimes(es)

get_dprimes <- function(es) {

    meta <- es$filt$mu
    names(meta) <- row.names(es$filt)

    contrasts <- list()
    dp_cols <- grep("^dp", colnames(es$raw), value = TRUE)

    for (col in dp_cols){
        dprimes <- es$raw[, col]
        names(dprimes) <- row.names(es$raw)
        contrasts[[col]] <- dprimes[!is.na(dprimes)]
    }
    return(list (meta=meta, contrasts=contrasts))
}


#---------------



# Divide integers up to n into n bins
#
# Function is useful to divide task among multiple threads.
#
# @param n number of items to divide into bins.
# @param nbins number of bins to divide n items into
#
# @return List of integer vectors specifying which items are in each bin
#
# @examples
# x <- 1:22268
# n <- 10
# bins <- split(x, sort(x%%n))

get_bins <- function(n, nbins) {
    split(1:n, sort(1:n%%nbins))
}


# Repeat vector in n rows of matrix.
#
# @param x Vector to repeat.
# @param n Number of rows to copy x into.
#
# @return Matrix with x in each of n rows.

rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
}
