# Setup drug combination training data
#
# User is asked to select drug1, drug2, and combo contrasts for each GSE.
# Selections are then used create training data.frame.
#
# @param diff_exprs Result from call to \code{\link[crossmeta]{diff_expr}}.
# @param prev_selections Used to re-use previous selections.
#
# @return data.frame

setup_combo_data <- function(diff_exprs, prev_selections=NULL) {

    #add moderated effect sizes
    diff_exprs <- add_dprime(diff_exprs)

    #get selections
    selections <- select_combo_data(diff_exprs, prev_selections)

    #combine
    data <- combine_combo_data(diff_exprs, selections)

    #return result
    combo_data <- list(data=data, selections=selections)
    return(combo_data)
}

#---------------

# Select drug1, drug2, and combo contrasts for each GSE
#
# Used by setup_combo_data to ask user to select drug1, drug2, and combo
# contrasts for each GSE.
#
# @inheritParams setup_combo_data
#
# @return list (one per GSE) of lists (one per combination treatment)
#    with contrast names for each combination treatment in each GSE.

select_combo_data <- function(diff_exprs, prev_selections){

    #setup selections
    selections <- vector("list", length(diff_exprs))
    names(selections) <- names(diff_exprs)

    #transfer prev_selections from GSEs also in diff_exprs
    prev_selections <- prev_selections[names(prev_selections) %in% names(selections)]
    selections[names(prev_selections)] <- prev_selections

    #select drug 1, drug 2, and combo from each study
    for (i in seq_along(diff_exprs)) {

        #check if previously selected
        gse_name <- names(diff_exprs)[i]
        if (!is.null(selections[[gse_name]])) next

        #if not, input selection
        choices <- names(diff_exprs[[i]]$top_tables)

        while (TRUE) {
            #select drug 1, 2, and combo contrast
            drug1 <- tcltk::tk_select.list(choices,
                                           title="Contrast for drug 1")
            if (drug1 == "") {break}
            drug2 <- tcltk::tk_select.list(choices,
                                           title="Contrast for drug 2")
            combo <- tcltk::tk_select.list(choices,
                                           title="Contrast for drug combo")

            #store selections
            selection <- c(drug1=drug1, drug2=drug2, combo=combo)
            selections[[gse_name]] <- append(selections[[gse_name]],
                                             list(selection))
        }
    }
    return(selections)
}

#---------------

# Creates final drug combination training data.frame
#
# Used by setup_combo_data to create final drug combination data.frame.
#
# @inheritParams setup_combo_data
#
# @return data.frame

combine_combo_data <- function(diff_exprs, selections){

    #return value
    data <- list()

    #loop through each GSE
    for (i in seq_along(selections)) {

        print(i)

        tops <- diff_exprs[[i]]$top_tables
        sels <- selections[[i]]
        ord <- row.names(tops[[1]])

        #loop through each selection within each GSE
        for (sel in sels) {

            drug1 <- tops[[ sel['drug1'] ]][ord, ]
            drug2 <- tops[[ sel['drug2'] ]][ord, ]
            combo <- tops[[ sel['combo'] ]][ord, ]

            #name and bind together
            colnames(drug1) <- paste("drug1", colnames(drug1), sep="_")
            colnames(drug2) <- paste("drug2", colnames(drug2), sep="_")
            colnames(combo) <- paste("combo", colnames(combo), sep="_")

            data[[length(data)+1]] <- cbind(drug1, drug2, combo, row.names=NULL)

            print(sum(is.na(data[[length(data)]]$combo_dprime)))
        }
    }
    data <- data.table::rbindlist(data)
    data <- as.data.frame(data)
    return(data)
}


#-----------------

# Add dprimes to top tables.
#
# Used by setup_combo_data to add moderated unbiased standardised effect sizes
# (dprimes) to top tables from differential expression analysis.
#
# @inheritParams setup_combo_data
#
# @return diff_exprs with dprimes added to top_tables for each contrast.

add_dprime <- function(diff_exprs) {

    for (study in names(diff_exprs)) {
        #get study degrees of freedom
        diff <- diff_exprs[[study]]
        df <- diff$ebayes_sv$df.residual + diff$ebayes_sv$df.prior

        for (con in names(diff$top_tables)) {
            #get sample sizes and top table for contrast
            classes <- pData(diff$eset)$treatment
            ni <- length(classes[classes == "ctrl"])
            nj <- length(classes[classes == "test"])

            tt <- diff$top_tables[[con]]

            #get dprime
            tt$dprime <- metaMA::effectsize(tt$t,
                                            ((ni * nj)/(ni + nj)),
                                            df)[, "dprime"]

            #store result
            diff_exprs[[study]]$top_tables[[con]] <- tt
        }
    }
    return(diff_exprs)
}


#-----------------



#' Make drug combination database.
#'
#' Creates an SQLite database with predicted effect sizes for all unique two-drug
#' combinations of drugs in the Connectivity Map build 2.
#'
#' Prediction model is a gradient boosted random forest trained on GEO microarray
#' studies where single treatments and their combinations were assayed.
#'
#' @import data.table ccdata
#' @param db_dir String specifying directory in which to create database.
#'    Default is working directory.
#'
#' @export
#' @seealso \link{query_combos}, \link{range_query_combos}.
#'
#' @examples \dontrun{
#' # make drug combination data base
#' make_drug_combos()
#' }

make_drug_combos <- function(db_dir=getwd()) {


    # Setup -------------------------


    #load model & cmap tables
    combo_model  <- ccdata::combo_model

    cmap_tables1 <- ccdata::cmap_tables1
    cmap_tables2 <- ccdata::cmap_tables2
    cmap_tables  <- c(cmap_tables1, cmap_tables2)

    drugs  <- names(cmap_tables)
    probes <- row.names(cmap_tables[[1]])

    #list of unique drug combos
    pairs <- utils::combn(drugs, 2, simplify=FALSE)

    #load annotation information
    hgu133a <- crossmeta:::get_biocpack("hgu133a.db")
    suppressMessages(map <- AnnotationDbi::select(hgu133a, probes, "SYMBOL"))
    map <- map[!is.na(map$SYMBOL),]
    map$SYMBOL <- toupper(map$SYMBOL)
    genes <- unique(map$SYMBOL)

    #make blank table in SQLite database
    db_path <- paste(db_dir, "drug_combos.sqlite", sep="/")
    db.pdc <- RSQLite::dbConnect(RSQLite::SQLite(), dbname=db_path)
    cols <- paste(shQuote(genes), "FLOAT", collapse=", ")

    statement <- paste("CREATE TABLE combo_preds",
                       "(drug_combo TEXT, ", cols, ")", sep="")

    RSQLite::dbSendQuery(db.pdc, statement)


    d1_cols <- paste("drug1", colnames(cmap_tables[[1]]), sep="_")
    d2_cols <- paste("drug2", colnames(cmap_tables[[1]]), sep="_")

    i <- 1
    pb  <- utils::txtProgressBar(min=1, max=length(pairs), style=3)
    for (pair in pairs){



        # Probe Predictions -------------------------


        dr1 <- pair[1]
        dr2 <- pair[2]

        #combine data from cmap_tables
        X <- cbind(cmap_tables[[dr1]], cmap_tables[[dr2]])
        X <- as.matrix(X)
        colnames(X) <- c(d1_cols, d2_cols)

        #make predictions with combo_model
        probe_preds <- xgboost::predict(combo_model, X)

        #predict using reverse order
        colnames(X) <- c(d2_cols, d1_cols)
        X <- X[, c(d1_cols, d2_cols)]
        probe_preds_rev <- xgboost::predict(combo_model, X)

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



        # Add to SQL Database -------------------------


        #format table
        gene_preds <- as.data.frame(gene_preds)
        row.names(gene_preds) <- gene_preds$SYMBOL
        gene_preds <- gene_preds[, -1, drop=FALSE]

        gene_preds <- as.data.frame(t(gene_preds))
        gene_preds$drug_combo <- row.names(gene_preds)
        row.names(gene_preds) <- NULL

        RSQLite::dbWriteTable(db.pdc, "combo_preds",
                              gene_preds[, c("drug_combo", genes)],
                              append = TRUE)

        utils::setTxtProgressBar(pb, i)
        i <- i + 1
    }
    close(pb)
}
