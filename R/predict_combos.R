#' Get predicted drug combination signatures.
#'
#' Predicts drug combination signatures.
#'
#' Drug combinations may more closely mimic or reverse the query signature than
#' individual drugs. The drug combinations that most closely mimic or reverse
#' a query signature can usually be determined by querying against drug
#' combination signatures predicted with the top few single drugs.
#'
#' @import data.table
#' @param with Character vector of drug names to combine cmap drugs with.
#'
#' @seealso \code{\link{query_drugs}} to determine overlap between query and
#'    predicted drug combination signatures.
#'
#' @return drug_info Matrix of differential expression values for drug
#'   combinations. Rows are genes, columns are drugs.
#' @export
#'
#' @examples
#' # this returns NULL:
#' combos_info <- predict_combos("cmap_drug")
#'
#'
#' # predicted expression values for combined treatment with sirolimus
#' # and all other cmap drugs:
#' # combos_info <- predict_combos("sirolimus")

predict_combos <- function(with) {

    # Setup -------------------------
    drugs <- colnames(ccdata::cmap_es)
    if (FALSE %in% (with %in% drugs))
        message("Drugs in 'with' not found."); return(NULL)

    #load model & cmap tables
    combo_model  <- ccdata::combo_model

    cmap_tables1 <- ccdata::cmap_tables1
    cmap_tables2 <- ccdata::cmap_tables2
    cmap_tables  <- c(cmap_tables1, cmap_tables2)

    drugs  <- drugs[!drugs %in% with]
    probes <- row.names(cmap_tables[[1]])

    #list of unique drug combos
    pairs <- expand.grid(with, drugs, stringsAsFactors = FALSE)

    #load annotation information
    hgu133a <- crossmeta:::get_biocpack("hgu133a.db")
    suppressMessages(map <- AnnotationDbi::select(hgu133a, probes, "SYMBOL"))
    map <- map[!is.na(map$SYMBOL),]
    map$SYMBOL <- toupper(map$SYMBOL)
    genes <- unique(map$SYMBOL)

    d1_cols <- paste("drug1", colnames(cmap_tables[[1]]), sep="_")
    d2_cols <- paste("drug2", colnames(cmap_tables[[1]]), sep="_")

    #final drug info data.frame
    drug_info <- data.frame(row.names = genes)

    i <- 1
    pb  <- utils::txtProgressBar(min=1, max=nrow(pairs), style=3)
    for (i in 1:nrow(pairs)) {

        # Probe Predictions -------------------------

        dr1 <- pairs[i, 1]
        dr2 <- pairs[i, 2]

        #combine data from cmap_tables
        X <- cbind(cmap_tables[[dr1]], cmap_tables[[dr2]])
        X <- as.matrix(X)
        colnames(X) <- c(d1_cols, d2_cols)

        #make predictions with combo_model
        probe_preds <- xgboost:::predict.xgb.Booster(combo_model, X)

        #predict using reverse order
        colnames(X) <- c(d2_cols, d1_cols)
        X <- X[, c(d1_cols, d2_cols)]
        probe_preds_rev <- xgboost:::predict.xgb.Booster(combo_model, X)

        #get average then multiply by 1.5 (model correction factor)
        probe_preds <- (probe_preds + probe_preds_rev) / 2 * 1.5

        #probe_preds vector to df
        probe_preds <- as.data.frame(probe_preds)
        colnames(probe_preds)  <- paste(dr1, dr2, sep=" + ")
        row.names(probe_preds) <- probes


        # Probes to Genes -------------------------


        #annotate probe probe_preds with SYMBOL
        probe_preds <- probe_preds[map$PROBEID, , drop=FALSE]
        probe_preds$SYMBOL <- map$SYMBOL

        #where duplicated SYMBOL, choose probe with largest absolute prediction
        probe_preds <- data.table(probe_preds)
        gene_preds  <- probe_preds[,
                                   lapply(.SD, function (col) col[which.max(abs(col))]),
                                   by='SYMBOL']


        # Add to Drug Info -------------------------

        gene_preds <- as.data.frame(gene_preds)
        row.names(gene_preds) <- gene_preds$SYMBOL
        gene_preds <- gene_preds[genes, -1, drop=FALSE]

        drug_info <- cbind(drug_info, gene_preds)
        utils::setTxtProgressBar(pb, i)
        i <- i + 1
    }
    return(as.matrix(drug_info))
}
