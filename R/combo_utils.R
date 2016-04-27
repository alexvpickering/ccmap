#' Title
#'
#' @import Biobase
#'
#' @param diff_exprs
#' @param prev_selections
#'
#' @return
#' @export
#'
#' @examples
setup_combo_data <- function(diff_exprs, prev_selections=NULL) {

    #add moderated effect sizes
    diff_exprs <- add_dprime(diff_exprs)

    #get selections
    selections <- select_combo_data(diff_exprs, prev_selections)

    #combine
    data <- combine_combo_data(diff_exprs, selections)

    #save probe order
    order <- featureNames(diff_exprs[[1]]$eset)

    #return result
    combo_data <- list(data=data, selections=selections, order=order)
    return(combo_data)
}

#---------------

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

combine_combo_data <- function(diff_exprs, selections){

    #return value
    data <- data.frame()
    ord <- row.names(diff_exprs[[1]]$eset)

    #loop through each GSE
    for (i in seq_along(selections)) {

        tops <- diff_exprs[[i]]$top_tables
        sels <- selections[[i]]

        #loop through each selection within each GSE
        for (sel in sels) {

            drug1 <- tops[[ sel['drug1'] ]][ord, ]
            drug2 <- tops[[ sel['drug2'] ]][ord, ]
            combo <- tops[[ sel['combo'] ]][ord, ]

            #name and bind together
            colnames(drug1) <- paste("drug1", colnames(drug1), sep="_")
            colnames(drug2) <- paste("drug2", colnames(drug2), sep="_")
            colnames(combo) <- paste("combo", colnames(combo), sep="_")

            data <- rbind(data, cbind(drug1, drug2, combo))
        }
    }
    return(data)
}


#-----------------

add_dprime <- function(diff_exprs) {

    for (study in names(diff_exprs)) {
        #get study degrees of freedom
        diff <- diff_exprs[[study]]
        df <- diff$ebayes_sv$df.residual + diff$ebayes_sv$df.prior

        for (con in names(diff$top_tables)) {
            #get sample sizes and top table for contrast
            classes <- diff$mama_data$clinicals[[con]]$treatment
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



#' Title
#'
#' @import data.table ccdata
#' @param db_dir
#'
#' @return
#' @export
#'
#' @examples
make_drug_combos <- function(db_dir=getwd()) {


    # Setup -------------------------


    #load model & cmap tables
    data(combo_model, package="ccdata")
    data(cmap_tables, package="ccdata")

    drugs  <- names(cmap_tables)
    probes <- row.names(cmap_tables[[1]])

    #list of unique drug combos
    pairs <- combn(drugs, 2, simplify=FALSE)

    #load annotation information
    suppressMessages(library("hgu133a.db"))
    suppressMessages(map <- AnnotationDbi::select(hgu133a.db, probes, "SYMBOL"))
    map <- map[!is.na(map$SYMBOL),]
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
    pb  <- txtProgressBar(min=1, max=length(pairs), style=3)
    for (pair in pairs){



        # Probe Predictions -------------------------


        dr1 <- pair[1]
        dr2 <- pair[2]

        #combine data from cmap_tables
        X <- cbind(cmap_tables[[dr1]], cmap_tables[[dr2]])
        X <- as.matrix(X)
        colnames(X) <- c(d1_cols, d2_cols)

        #make predictions with combo_model
        probe_preds <- xgboost::predict(combo_model, X) - 0.5

        #predict using reverse order
        colnames(X) <- c(d2_cols, d1_cols)
        X <- X[, c(d1_cols, d2_cols)]
        probe_preds_rev <- xgboost::predict(combo_model, X) - 0.5

        #get average
        probe_preds <- (probe_preds + probe_preds_rev) / 2

        #probe_preds vector to df
        probe_preds <- as.data.frame(probe_preds)
        colnames(probe_preds)  <- paste(dr1, dr2, sep=" + ")
        row.names(probe_preds) <- probes



        # Probes to Genes -------------------------


        #annotate probe probe_preds with SYMBOL
        probe_preds <- probe_preds[map$PROBEID, , drop=FALSE]
        probe_preds$SYMBOL <- map$SYMBOL

        #where duplicated SYMBOL, choose probe with prediction furthest from 0.5
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

        setTxtProgressBar(pb, i)
        i <- i + 1
    }
    close(pb)
}
