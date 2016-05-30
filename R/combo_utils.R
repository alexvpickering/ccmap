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
    diff_exprs <- add_es(diff_exprs, cols = "dprime")

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
        }
    }
    data <- data.table::rbindlist(data)
    data <- as.data.frame(data)
    return(data)
}


#-----------------

# Add metaMA effectsize values to top tables.
#
# Used by setup_combo_data to add moderated unbiased standardised effect sizes
# (dprimes) to top tables from differential expression analysis.
#
# @inheritParams setup_combo_data
#
# @return diff_exprs with dprimes added to top_tables for each contrast.


add_es <- function(diff_exprs, cols = c("dprime", "vardprime")) {

    for (i in seq_along(diff_exprs)) {

        # get study degrees of freedom and group classes
        study <- diff_exprs[[i]]

        df <- study$ebayes_sv$df.residual + study$ebayes_sv$df.prior
        classes <- Biobase::pData(study$eset)$group

        for (con in names(study$top_tables)) {
            # get group names for contrast
            groups <- gsub("GSE.+?_", "", con)
            groups <- strsplit(groups, "-")[[1]]

            # get sample sizes for groups
            ni <- sum(classes == groups[2])
            nj <- sum(classes == groups[1])

            # bind effect size values with top table
            tt <- study$top_tables[[con]]
            es <- metaMA::effectsize(tt$t, ((ni * nj)/(ni + nj)), df)[, cols, drop = FALSE]
            tt <- cbind(tt, es)

            # store results
            study$top_tables[[con]] <- tt
        }
        diff_exprs[[i]] <- study
    }
    return(diff_exprs)
}
