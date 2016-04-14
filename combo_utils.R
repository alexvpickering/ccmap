library(metaMA)
library(reshape2)

setup_combo_data <- function(diff_exprs, prev_selections=NULL, by="array") {

    #add moderated effect sizes
    diff_exprs <- add_dprime(diff_exprs)

    #get selections
    selections <- select_combo_data(diff_exprs, prev_selections)

    #combine
    data <- combine_combo_data(diff_exprs, selections, by)

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
            drug1 <- tk_select.list(choices, title="Contrast for drug 1")
            if (drug1 == "") {break}
            drug2 <- tk_select.list(choices, title="Contrast for drug 2")
            combo <- tk_select.list(choices, title="Contrast for drug combo")

            #store selections
            selection <- c(drug1=drug1, drug2=drug2, combo=combo)
            selections[[gse_name]] <- append(selections[[gse_name]], list(selection))
        }
    }
    return(selections)
}

#---------------

combine_combo_data <- function(diff_exprs, selections, by="array"){

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

            if (by == "array") data <- by_array(drug1, drug2, combo, data)
            if (by == "probe") data <- rbind(data, cbind(drug1, drug2, combo))

        }
    }
    return(data)
}

by_array <- function(drug1, drug2, combo=NULL, data=data.frame()) {
    #OUT: samples (rows) all probes for drug1, drug2, and combo

    #add probes column
    drug1$probe <- row.names(drug1)
    drug2$probe <- row.names(drug2)

    #cast to single row
    drug1 <- dcast(melt(drug1, id.var="probe"), 1 ~ variable + probe)[, -1]
    drug2 <- dcast(melt(drug2, id.var="probe"), 1 ~ variable + probe)[, -1]


    if (!is.null(combo)) {

        #same for combo data
        combo$probe <- row.names(combo)
        combo <- dcast(melt(combo, id.var="probe"), 1 ~ variable + probe)[, -1]

        #put together and add to data
        data <- rbind(data, cbind(drug1, drug2, combo))

    } else data <- rbind(data, cbind(drug1, drug2))
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
            tt$dprime <- effectsize(tt$t, ((ni * nj)/(ni + nj)), df)[, "dprime"]

            #store result
            diff_exprs[[study]]$top_tables[[con]] <- tt
        }
    }
    return(diff_exprs)
}