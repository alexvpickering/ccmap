
#' Title
#'
#' @importFrom foreach foreach %dopar%
#' @param query_genes
#' @param db_dir
#' @param query_n
#' @param drug_n
#'
#' @return
#' @export
#'
#' @examples
query_combos <- function(query_genes, db_dir,
                        query_n=length(query_genes),
                        ncores=parallel::detectCores(),
                        step=NULL) {

    # 856086 combos = 654 * 1309
    doMC::registerDoMC(ncores)

    suppressWarnings(bins <- split(1:654, 1:ncores))
    bins <- sapply(bins, length)

    res <- foreach(i=1:ncores) %dopar% {

        #connect to db
        db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname=db_dir)

        #a <- sum(bins[1:i]) - bins[i] + 1
        #b <- sum(bins[1:i])

        a <- i*2-1
        b <- i*2

        tmp <- list()
        for (j in a:b) {

            #get preds for drug combos
            statement <- paste("SELECT * from combo_preds WHERE rowid BETWEEN",
                               (j*1309)-1308, "AND", j*1309)

            combo_preds <- RSQLite::dbGetQuery(db, statement)

            #format predictions
            combo_names <- combo_preds$drug_combo
            combo_preds <- as.data.frame(t(combo_preds[,-1]))
            colnames(combo_preds)  <- combo_names

            #query combo predictions
            combos_res <- query_drugs(query_genes, as.matrix(combo_preds),
                                     query_n, step=step)

            #store result
            tmp[[ length(tmp)+1 ]] <- combos_res
        }
        RSQLite::dbDisconnect(db)
        unlist(tmp, recursive=FALSE)
    }

    #gather results with same query size
    res <- unlist(res, recursive=FALSE)
    ns  <- unique(names(res))
    res <- lapply(ns, function(n) clean_rbindlist(res[names(res) == n]))
    names(res) <- ns

    #sort results by net overlap
    res <- lapply(res, function(x) x[order(x$net, decreasing=T), ])
    return(res)
}


#---------------------


#' Title
#'
#' @param query_genes
#' @param db_dir
#' @param ncores
#'
#' @return
#' @export
#'
#' @examples

range_query_combos <- function(query_genes, db_dir,
                              ncores=parallel::detectCores(), step=100) {

    #query combos for a range of query sizes
    combos_res <- query_combos(query_genes, db_dir, ncores=ncores, step=step)

    #create ranks and nets dataframes
    return(get_range_res(combos_res))
}
