#' Get overlap between query and predicted drug combination signatures.
#'
#' Drugs with the largest positive and negative cosine similarity are predicted to,
#' respectively, mimic and reverse the query signature. Values range from +1 to -1.
#'
#' To predict and query all 856086 two-drug cmap combinations, the 'average'
#' \code{method} can take as little as 10 minutes (Intel Core i7-6700). The 'ml'
#' (machine learning) \code{method} takes two hours on the same hardware and
#' requires ~10GB of RAM but is slightly more accurate. Both methods will run
#' faster by specifying only a subset of drugs using the \code{include} parameter.
#' To speed up the 'ml' method, the MRO+MKL distribution of R can help
#' substantially (\href{https://mran.revolutionanalytics.com/open/}{link}).
#' The combinations of LINCS l1000 signatures (~26 billion) can also be queried
#' using the 'average' \code{method}. In order to compare l1000 results to those
#' obtained with cmap, only the same genes should be queried (see example).
#'
#' @importFrom foreach foreach %dopar%
#'
#' @param query_genes Named numeric vector of differentual expression values for
#'   query genes. Usually 'meta' slot of \code{get_dprimes} result.
#' @param drug_info Character vector specifying which dataset to query
#'   (either 'cmap' or 'l1000'). Can also provide a matrix of differential expression
#'   values for drugs or drug combinations (rows are genes, columns are drugs).
#' @param method One of 'average' (default) or 'ml' (machine learning -
#'   see details and vignette).
#' @param include Character vector of drug names for which combinations with all
#'   other drugs will be predicted and queried. If \code{NULL} (default),
#'   all two drug combinations will be predicted and queried.
#' @param ncores Integer, number of cores to use for method 'average'. Default is
#'   to use all cores.
#'
#' @return Vector of cosine similarities between query and drug combination signatures.
#' @export
#'
#' @examples
#' library(lydata)
#' library(crossmeta)
#'
#' # location of data
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' # gather GSE names
#' gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")
#'
#' # load previous analysis
#' anals <- load_diff(gse_names, data_dir)
#'
#' # perform meta-analysis
#' es <- es_meta(anals)
#'
#' # get dprimes
#' dprimes <- get_dprimes(es)
#'
#' # query combinations of metformin and all other cmap drugs
#' top_met_combos <- query_combos(dprimes$meta, include = 'metformin', ncores = 1)
#'
#' # previous query but with machine learning method
#' # top_met_combos <- query_combos(dprimes$meta, method = 'ml', include = 'metformin')
#'
#' # query all cmap drug combinations
#' # top_combos <- query_combos(dprimes$meta)
#'
#' # query all cmap drug combinations with machine learning method
#' # top_combos <- query_combos(dprimes$meta, method = 'ml')
#'
#' # query l1000 and cmap using same genes
#' # library(ccdata)
#' # data(cmap_es)
#' # data(l1000_es)
#' # cmap_es <- cmap_es[row.names(l1000_es), ]
#'
#' # met_cmap  <- query_combos(dprimes$meta, cmap_es,  include = 'metformin')
#' # met_l1000 <- query_combos(dprimes$meta, l1000_es, include = 'metformin')

query_combos <- function(query_genes, drug_info = c('cmap', 'l1000'), method = c('average', 'ml'), include = NULL, ncores=parallel::detectCores()) {

    if (method[1] == 'ml' & drug_info[1] != 'cmap') {
        stop("Machine learning method only available for 'cmap' dataset.")
    }

    if (class(drug_info) == 'character') {
        # default to cmap_es for drug_info
        fname <- paste0(drug_info[1], '_es')
        utils::data(list = fname, package = "ccdata", envir = environment())
        drug_info <- get(fname)
        rm(list = fname)
    }

    drugs   <- colnames(drug_info)
    cpds    <- gsub('_.+?$', '', drugs)
    include <- drugs[cpds %in% include]

    # check 'include'
    if (FALSE %in% (include %in% drugs)) {
        message("Drugs in 'include' not found.")
        return(NULL)
    }

    # query all combinations if include is NULL
    if (is.null(include)) include <- drugs

    # use average model
    if (method[1] != "ml") {
        return(query_combos_average(query_genes, drug_info, include, ncores))

    # use machine learning 'ml' model
    } else {
        dat  <- load_ccdata(drug_info)

        exclude <- c()
        res <- c()
        i   <- 1
        pb  <- utils::txtProgressBar(min=1, max=length(include), style=3)

        for (drug in include) {
            combos_es <- predict_combos(drug, exclude, dat)
            if (is.null(combos_es)) break

            res <- c(res, query_drugs(query_genes, combos_es, sorted = FALSE))

            exclude <- c(exclude, drug)
            utils::setTxtProgressBar(pb, i)
            i <- i + 1
        }
        return(sort(res, decreasing = TRUE))
    }
}



#---------------



# Query drug combinations using average model.
#
# @param query_genes Named numeric vector of differentual expression values for
#   query genes. Usually 'meta' slot of \code{get_dprimes} result.
# @param cmap_es Matrix with cmap dprime values.
# @param include Character vector of cmap drug names for which combinations with all
#   other cmap drugs will be predicted and queried. If \code{NULL} (default),
#   all 856086 two drug combinations will be predicted and queried.
#
# @return Vector of numeric values between 1 and -1 indicating extent of overlap
#   between query and drug combination signatures.

query_combos_average <- function(query_genes, drug_info, include, ncores) {
    i = NULL  # bind global variable

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    drugs <- colnames(drug_info)
    bins  <- get_bins(length(include), ncores)

    resl <- foreach::foreach(i=1:min(ncores, length(bins))) %dopar% {

        bin <- bins[[i]]

        #add to exclude list (for all but first bin)
        if (i != 1) {
            exclude <- drugs[1:(bin[1]-1)]
        } else {
            exclude <- c()
        }

        res <- c()
        for (drug in include[bin]) {

            # average drug and all other drugs (except excluded)
            other_drugs <- setdiff(drugs, c(exclude, drug))
            if (length(other_drugs) == 0) break
            drug_es <- drug_info[, drug]

            # break up task up to reduce memory footprint
            dbins <- get_bins(length(other_drugs), 10)

            for (dbin in dbins) {
                combo_es <- drug_info[, other_drugs[dbin], drop = FALSE]
                combo_es <- sweep(combo_es, 1, drug_es, `+`)/2

                # update colnames to reflect combo
                colnames(combo_es) <- paste(drug, colnames(combo_es), sep = " + ")

                # query combo_es
                res <- c(res, query_drugs(query_genes, combo_es, FALSE))
            }
            # add drug to excluded list
            exclude <- c(exclude, drug)
        }
        res
    }
    parallel::stopCluster(cl)
    return(sort(unlist(resl), decreasing = TRUE))
}




#---------------




# Load data for machine learning method of predict_combos.
#
# @param cmap_es Matrix with cmap dprime values. If NULL, will be loaded.
#
# @return List with objects need for 'ml' method of predict_combos

load_ccdata <- function(cmap_es = NULL) {

    #bind global variables
    cmap_var = net1 = net2 = genes = xgb_mod = NULL

    if(is.null(cmap_es))
        utils::data("cmap_es", package = "ccdata", envir = environment())

    # load everything
    utils::data("cmap_var", "net1", "net2", "genes", "xgb_mod",
                package = "ccdata", envir = environment())

    # keep/order genes as for nnet
    cmap_es  <- t(cmap_es[genes, ])
    cmap_var <-  cmap_var[genes, ]

    return(list(es=cmap_es, var=cmap_var, net1=net1, net2=net2,
                genes=genes, xgb_mod=xgb_mod))
}
