#' Title
#'
#' @importFrom foreach foreach %dopar%
#'
#' @param ncores
#' @param query_genes
#'
#' @return
#' @export
#'
#' @examples

query_combos <- function(query_genes, method = "avg", include = NULL) {

    # get cmap_es
    utils::data("cmap_es", package = "ccdata", envir = environment())
    drugs <- colnames(cmap_es)

    # check 'include'
    if (FALSE %in% (include %in% drugs)) {
        message("Drugs in 'include' not found.")
        return(NULL)
    }

    # query all combinations if include is NULL
    if (is.null(include)) include <- drugs

    # use average model
    if (method != "ml") {
        return(query_combos_avg(query_genes, cmap_es, include))

    # use machine learning 'ml' model
    } else {
        dat  <- load_data(cmap_es)

        exclude <- c()
        res     <- c()
        for (drug in include) {
            combos_es <- predict_combos(drug, exclude, dat)
            res <- c(res, query_drugs(query_genes, combos_es, sorted = FALSE))
            exclude <- c(exclude, drug)
        }
        return(sort(res, decreasing = TRUE))
    }
}




query_combos_avg <- function(query_genes, cmap_es, include) {

    ncores=parallel::detectCores()
    doMC::registerDoMC(ncores)

    drugs <- clnames(cmap_es)
    bins  <- get_bins(length(include), ncores)

    resl <- foreach(i=1:min(ncores, length(bins))) %dopar% {

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
            drug_es <- cmap_es[, drug]

            # break up task up to reduce memory footprint
            dbins <- avpick::get_bins(length(other_drugs), 10)

            for (dbin in dbins) {
                combo_es <- cmap_es[, other_drugs[dbin], drop = FALSE]
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
    return(sort(unlist(resl), decreasing = TRUE))
}



load_data <- function(cmap_es = NULL) {

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
