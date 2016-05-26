#' Get overlap between query and drug combination signatures.
#'
#' Determines the number of genes that change in the same direction between
#' query and predicted drug combination signatures.
#'
#' All 856086 combinations of the 1309 drugs from the Connectivity Map build 02
#' release are compared to the query signature. Drug combinations with the
#' largest positive and negative net overlap are predicted to, respectively,
#' mimic and reverse the query signature. Drug combinations may more closely
#' mimic or reverse the query signature than individual drugs.
#'
#' @importFrom foreach foreach %dopar%
#'
#' @inheritParams query_drugs
#' @param db_dir String specifying full path to drug_combos.sqlite database.
#' @param ncores Number of cores to use for parallel queries. Default is all
#'    available.
#'
#' @family query functions
#' @seealso \code{\link{make_drug_combos}} to create drug combination database.
#'
#' @return List of data.frames (one per query size) each with columns:
#'   \item{drug_n}{Number of drug combination genes used.}
#'   \item{query_n}{Number of query genes used.}
#'   \item{overlap}{Number of genes that change in the same direction in drug
#'      combination and query signatures.}
#'   \item{cross}{Number of genes that change in the opposite direction in drug
#'      combination and query signatures.}
#'   \item{net}{Difference between overlap and cross.}
#' @export
#'
#' @examples \dontrun{
#' library(crossmeta)
#' library(lydata)
#'
#' db_dir <- "path/to/drug_combos.sqlite"
#'
#' # location of data
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' # gather GSE names
#' gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")
#'
#' # load previous analysis
#' anals <- load_diff(gse_names, data_dir)
#'
#' # perform meta-analysis
#' es <- es_meta(anals)
#'
#' # get query signature
#' dprimes <- get_dprimes(es)
#'
#' # query drug combination database
#' top_combos <- query_combos(dprimes$meta, db_dir)
#' }

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

        a <- sum(bins[1:i]) - bins[i] + 1
        b <- sum(bins[1:i])

        #a <- i*2-1
        #b <- i*2

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
    res <- lapply(res, function(x) x[order(x$net, decreasing = TRUE), ])
    return(res)
}


#---------------------


#' Get overlap between query and drug combination signatures for a range of
#' query gene sizes.
#'
#' Queries drug combinations database using a range of query gene sizes. Results
#' are sorted by an auc metric which has the advantage of weighting query genes
#' according to their extent of differential expression.
#'
#' All 856086 combinations of the 1309 drugs from the Connectivity Map build 02
#' release are compared to the query signature. Drug combinations with the
#' largest positive and negative net overlap are predicted to, respectively,
#' mimic and reverse the query signature. Drug combinations may more closely
#' mimic or reverse the query signature than individual drugs.
#'
#' @inheritParams query_combos
#' @family query functions
#'
#' @return data.frame with number of net genes that overlap between query and drug
#'    combination signatures. Columns correspond to query sizes, rows to drug
#'    combinations. Combinations are sorted by decreasing area under the curve
#'    formed by plotting net overlap as a function of query size.
#' @export
#'
#' @examples \dontrun{
#' library(crossmeta)
#' library(lydata)
#'
#' db_dir <- "path/to/drug_combos.sqlite"
#'
#' # location of data
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' # gather GSE names
#' gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")
#'
#' # load previous analysis
#' anals <- load_diff(gse_names, data_dir)
#'
#' # perform meta-analysis
#' es <- es_meta(anals)
#'
#' # get query signature
#' dprimes <- get_dprimes(es)
#'
#' # query drug combination database for range of query sizes
#' ranges_res <- range_query_combos(dprimes$meta, db_dir)
#' }

range_query_combos <- function(query_genes, db_dir,
                               ncores=parallel::detectCores(), step=100) {

    #query combos for a range of query sizes
    combos_res <- query_combos(query_genes, db_dir, ncores=ncores, step=step)

    #create ranks and nets dataframes
    return(get_range_res(combos_res))
}
