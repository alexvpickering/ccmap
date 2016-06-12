# Test meta-analysis against contrasts for a drug in CMAP.
#
# Use for gold-standard evaluation of crossmeta. Namely, does a meta-analysis
# using independent microarray data for drugs present in CMAP improve query
# results as compared to the contrasts used in the meta-analysis.
#
# @param dprimes Result of call to \code{get_dprimes}.
#
# @return List with results for meta-analysis and contrasts.

test_ma <- function(dprimes, suffix = "") {

    # get contrast results for range of query sizes
    cons_res <- lapply(dprimes$contrasts, function(con) {
        query_drugs(con, step = 100)
    })

    # get meta-analysis results for range of query sizes
    meta_res <- query_drugs(dprimes$meta, step = 100)

    res <- list(meta = meta_res, cons = cons_res)
    saveRDS(res, paste0("ma_res", suffix, ".rds"))

    return(res)
}
