# Get predicted drug combination signatures.
#
# Predicts drug combination signatures.
#
# Drug combinations may more closely mimic or reverse the query signature than
# individual drugs. The drug combinations that most closely mimic or reverse
# a query signature can usually be determined by querying against drug
# combination signatures predicted with the top few single drugs.
#
# @import data.table ccdata AnnotationDbi BiocInstaller xgboost
# @param include Character vector of drug names to combine cmap drugs with.
# @param exclude Character vector of drug names to not combine with include.
# @param dat List with data for machine learning model. If NULL (default), will
#   be loaded automatically.
#
# @seealso \code{\link{query_drugs}} to determine overlap between query and
#    predicted drug combination signatures.
#
# @return drug_info Matrix of differential expression values for drug
#   combinations. Rows are genes, columns are drugs.
#
# @examples
# # predicted expression values for combined treatment with sirolimus
# # and all other cmap drugs:
# # combos_info <- predict_combos("sirolimus")

predict_combos <- function(include, exclude = NULL, dat = NULL) {

    if (is.null(dat)) {
        dat <- load_ccdata()
    }

    # unpack dat
    cmap_es  <- dat$es
    cmap_var <- dat$var
    net1     <- dat$net1
    net2     <- dat$net2
    genes    <- dat$genes
    xgb_mod  <- dat$xgb_mod
    rm(dat)

    res <- list()
    for (drug in include) {

        # predict combinations of 'include' and other drugs except 'exclude'
        other_drugs <- setdiff(row.names(cmap_es), c(drug, exclude))

        if (length(other_drugs) == 0) return(NULL)

        # setup test data for NNet predictions
        drug_es <- rep.row(cmap_es[drug, ], length(other_drugs))
        Xnet    <- cbind(drug_es, cmap_es[other_drugs, , drop=FALSE])
        rm(drug_es)

        # average NNet predictions from net1 and net2
        preds <- predict.nnet(net1, Xnet)
        preds <- (preds + predict.nnet(net2, Xnet))/2

        # setup test data for xgb predictions
        Xgb <- matrix(nrow = length(preds), ncol = 5)
        colnames(Xgb) <- c("net_preds", "drug1_dprime", "drug2_dprime",
                            "drug1_vardprime", "drug2_vardprime")

        # add NNet preds and dprimes
        Xgb[, 1] <- as.vector(t(preds))
        rm(preds)

        Xnet1 <- t(Xnet[, 1:11525])
        Xnet2 <- t(Xnet[, 11526:23050]) ; rm(Xnet)
        Xgb[, 2:3] <- cbind(as.vector(Xnet1), as.vector(Xnet2))
        rm(Xnet1, Xnet2)

        # add vardprimes
        Xgb[, 4] <- cmap_var[, drug]
        Xgb[, 5] <- as.vector(cmap_var[, other_drugs])

        # xgboost predictions
        combos_es <- xgboost:::predict.xgb.Booster(xgb_mod, Xgb)
        rm(Xgb)

        # setup combos_es
        dim(combos_es) <- c(11525, length(combos_es)/11525)
        colnames(combos_es)  <- paste(drug, other_drugs, sep = " + ")
        row.names(combos_es) <- genes

        res[[length(res)+1]] <- combos_es
        rm(combos_es)

        exclude <- c(exclude, drug)
    }
    return(do.call(cbind, res))
}



#---------------






# Downloads bioconductor package.
#
# Used by symbol_annot to download annotation data packages from bioconductor.
#
# @param biocpack_name String specifying bioconductor package to download.
#
# @seealso \link{get_biocpack_name}, \link{symbol_annot}.
# @return NULL (downloads and loads requested package).

get_biocpack <- function(biocpack_name) {

    if (!requireNamespace(biocpack_name, quietly = TRUE)) {
        BiocInstaller::biocLite(biocpack_name)
    }
    db <- get(biocpack_name, getNamespace(biocpack_name))
    return (db)
}



#---------------




# Get neural network predictions.
#
# Neural network has one hidden layer and very_leaky_rectifier activation
# (hence pmax(z2/3, z2)).
#
# @param net List with W1/W2 weight matrices and b1/b2 bias vectors.
# @param X Matrix to get predictions for.
#
predict.nnet <- function(net, X) {
    z2 <- sweep(X %*% net$W1, 2, net$b1, "+")
    a2 <- pmax(z2/3, z2)
    return (sweep(a2 %*% net$W2, 2, net$b2, "+"))
}





