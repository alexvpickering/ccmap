library(metaMA)
library(reshape2)

setup_combo_data <- function(diff_exprs, prev_selections=NULL, by="array") {

    #add moderated effect sizes
    diff_exprs <- add_dprime(diff_exprs)

    #get selections
    selections <- select_combo_data(diff_exprs, prev_selections)

    #combine
    data <- combine_combo_data(diff_exprs, selections, by)

    #return result
    combo_data <- list(data=data, selections=selections)
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

by_array <- function(drug1, drug2, combo, data){
    #OUT: samples (rows) all probes for drug1, drug2, and combo

    #add probes column
    drug1$probe <- row.names(drug1)
    drug2$probe <- row.names(drug2)
    combo$probe <- row.names(combo)

    #cast to single row
    drug1 <- dcast(melt(drug1, id.var="probe"), 1 ~ variable + probe)[, -1]
    drug2 <- dcast(melt(drug2, id.var="probe"), 1 ~ variable + probe)[, -1]
    combo <- dcast(melt(combo, id.var="probe"), 1 ~ variable + probe)[, -1]

    #put together and add to data
    data <- rbind(data, cbind(drug1, drug2, combo))

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


#------------------

# add_top_cors <- function(combo_probe, probes_cor, n=1) {
#   #adds info from top n most correlated probes

#   #only consider probes in both combo_probe and probes_cor
#   probes <- colnames(probes_cor)
#   probes <- probes[probes %in% combo_probe$keys]
#   probes_cor <- probes_cor[probes, probes]

#   #identify top n correlating probes for each probe
#   top_cors <- list()
  
#   for (probe in probes) {
#     probe_cor <- sort(probes_cor[, probe], decreasing=T)
#     top_cors[[probe]] <- probe_cor[2:(2+n-1)]
#   }

#   #add dprime, adj.P.Val, and correlation for drug1 and drug2
#   for (i in 1:n) {

#     #new column names(dprime, pval, and cor)
#     dp1_name <- paste("drug1_", "cor", i, "_dprime", sep="")
#     dp2_name <- paste("drug2_", "cor", i, "_dprime", sep="")

#     p1_name <- paste("drug1_", "cor", i, "_adj.P.Val", sep="")
#     p2_name <- paste("drug2_", "cor", i, "_adj.P.Val", sep="")

#     c_name <- paste("cor", i, sep="")

#     #prep
#     probe_rows <- paste(rep(probes, each=80), c("", 1:79), sep="")
#     dp1 <- c()
#     dp2 <- c()
#     p1  <- c()
#     p2  <- c()
#     c   <- c()

#     for (probe in probes){

#       #get probe (and corresponding rows) with top_i correlation
#       pci <- top_cors[[probe]][i]
#       pci_rows <- paste(names(pci), c("", 1:79), sep="")

#       #get dprime of pci for drug1 and drug2
#       dp1 <- c(dp1, combo_probe[pci_rows, drug1_dprime])
#       dp2 <- c(dp2, combo_probe[pci_rows, drug2_dprime])

#       #get adj.P.Val of pci for drug1 and drug2
#       p1 <- c(p1, combo_probe[pci_rows, drug1_adj.P.Val])
#       p2 <- c(p2, combo_probe[pci_rows, drug2_adj.P.Val])

#       #add correlation value of pci to rows for probe
#       c <- c(c, rep(pci, 80))
#     }

#     #add values to combo_probe
#     combo_probe[probe_rows, eval(dp1_name)] <- dp1
#     combo_probe[probe_rows, eval(dp2_name)] <- dp2

#     combo_probe[probe_rows, eval(p1_name)] <- p1
#     combo_probe[probe_rows, eval(p2_name)] <- p2

#     combo_probe[probe_rows, eval(c_name)] <- c
#   }
#   return(combo_probe)
# }