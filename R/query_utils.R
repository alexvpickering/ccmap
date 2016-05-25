#' Extract unbiased effect sizes from meta-analysis by crossmeta.
#'
#' Function extracts mu (overall mean effect size) and dprimes (unbiased
#' effect sizes from each contrast).
#'
#' Result used to query connectivity map drugs and predicted drug combinations.
#'
#' @param es Result of call to \code{es_meta}.
#'
#' @return List containing:
#'   \item{meta}{Named numeric vector with overall mean effect sizes for all genes
#'      from meta-analysis.}
#'   \item{contrasts}{List of named numeric vectors (one per contrast) with
#'      unbiased effect sizes for all measured genes.}
#' @export
#'
#' @seealso \code{\link{es_meta}}.
#' @examples
#' library(crossmeta)
#' library(lydata)
#'
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' # gather GSE names
#' gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")
#'
#' # load previous differential expression analysis
#' anals <- load_diff(gse_names, data_dir)
#'
#' # run meta-analysis
#' es <- es_meta(anals)
#'
#' #get dprimes
#' dprimes <- get_dprimes(es)

get_dprimes <- function(es) {

    meta <- es$filt$mu
    names(meta) <- row.names(es$filt)

    contrasts <- list()
    dp_cols <- grep("^dp", colnames(es$raw), value = TRUE)

    for (col in dp_cols){
        dprimes <- es$raw[, col]
        names(dprimes) <- row.names(es$raw)
        contrasts[[col]] <- dprimes[!is.na(dprimes)]
    }
    return(list (meta=meta, contrasts=contrasts))
}


#---------------------


# Get list of top differentially regulated genes for drugs.
#
# @inheritParams query_drugs
#
# @seealso query_drugs
# @return List of lists (one per drug in drug_info) with slots:
#   \item{up}{Named vector of genes up-regulated by drug.}
#   \item{up}{Named vector of genes down-regulated by drug.}

get_drug_genes <- function(drug_info, query_genes,
                           drug_n=nrow(drug_info), es=TRUE) {

    #fix for my drug_combos.sqlite (remove if re-create)
    row.names(drug_info) <- toupper(row.names(drug_info))

    drug_genes <- list()
    drugs <- colnames(drug_info)

    #only consider drug genes that are also in query genes
    drug_info <- drug_info[row.names(drug_info) %in%
                               toupper(names(query_genes)), ]

    #max drug_n is drug genes that are also in query genes
    if (drug_n > nrow(drug_info)) {
        drug_n <- nrow(drug_info)
    }
    for (drug in drugs) {
        if (es) {
            #get up/dn for each drug
            dr <- drug_info[, drug]
            updn <- get_updn(dr, drug_n)

            up <- updn$up
            dn <- updn$dn
        } else {
            dr <- sort(drug_info[, drug])

            up <- names (head(dr, drug_n / 2))
            dn <- names (tail(dr, drug_n / 2))
        }
        #up/dn regulated
        top_list <- list(up=up, dn=dn)
        drug_genes[[drug]] <- top_list
    }
    return(drug_genes)
}

#----------------

# Get top DE genes for query.
#
# @inheritParams query_drugs
#
# @return List with slots 'up' and 'dn' each containing captalized character
#   vectors of up- and down-regulated genes.

setup_query_genes <- function(query_genes, query_n, drug_info) {

    #make names uppercase
    names(query_genes) <- toupper(names(query_genes))

    #remove genes not in drug_info
    query_genes <- query_genes[names(query_genes) %in% row.names(drug_info)]

    updn <- get_updn(query_genes, query_n)

    query_list <- list(up=updn$up, dn=updn$dn)
    return(query_list)
}

#---------------------

# Get number of top DE genes.
#
# Function returns at most the requested number of genes, attempting to devide
# equally between the top up and down regulated genes.
#
# If fewer than half of requested genes are available among up (or down)
# regulated genes, then more genes will be taken from down (or up) regulated
# genes.
#
# @inheritParams query_drugs
#
# @seealso \code{\link{setup_query_genes}}, \code{\link{get_drug_genes}}
# @return List with slots 'up' and 'dn' each containing character vectors
#   of up- and down-regulated genes.

get_updn <- function (genes, genes_n) {

    #sort query genes
    genes <- sort(genes, decreasing=TRUE)

    #get up/dn list
    up <- names(genes[genes > 0])
    dn <- names(genes[genes < 0])

    if (length(up) < genes_n / 2 &
        length(dn) > genes_n / 2) {
        #take more genes from dn
        n_up <- length(up)
        n_dn <- genes_n - length(up)
    } else if (length(dn) < genes_n / 2 &
               length(up) > genes_n / 2) {
        #take more genes from up
        n_up <- genes_n - length(dn)
        n_dn <- length(dn)
    } else {
        #take genes_n / 2 from up and dn
        n_up <- genes_n / 2
        n_dn <- genes_n / 2
    }
    #get n_up, n_dn from up/dn list
    up <- head(up, n_up)
    dn <- tail(dn, n_dn)

    return(list(up=up, dn=dn))
}

#---------------------

# Get overlap between drug and query signatures.
#
# Determines number of genes regulated in the same direction between drug and
# query signatures.
#
# @inheritParams query_drugs
#
# @return Data.frame with query_n, drug_n, and overlap columns.

get_overlap <- function(drug_genes, query_genes) {
    #used to provide overlap between drug_genes and query_genes

    drugs <- names(drug_genes)
    overlap <- data.frame(row.names=drugs, stringsAsFactors=FALSE)

    #add number of drug/query genes
    overlap$query_n <- length(query_genes$up) + length(query_genes$dn)
    overlap$drug_n  <- length(drug_genes[[1]]$up) + length(drug_genes[[1]]$dn)

    for (dr in drugs) {
        dr_genes <- drug_genes[[dr]]

        #get up/down and total overlap
        up <- intersect(dr_genes$up, query_genes$up)
        dn <- intersect(dr_genes$dn, query_genes$dn)
        tot <- length(up) + length(dn)

        #add total overlap for drug
        overlap[dr, "overlap"] <- tot

    }
    return (overlap)
}



#---------------------



# Format results for range_query functions.
#
# Generates data.frame with net overlap at each query size for all drugs or drug
# combinations. Results are sorted by auc from nets vs query gene sizes.
#
# @param query_res Reslt of call to query_drugs

get_range_res <- function(query_res) {

    # bind 'net' from data.frames for all query sizes
    query_ns <- as.integer(names(query_res))
    drugs <- row.names(query_res[[1]])

    query_res <- lapply(query_res, function(x) x[drugs, "net", drop = FALSE])
    query_res <- do.call(cbind, query_res)
    colnames(query_res) <- query_ns

    #sort query_res by auc formed by plot of nets vs query_n
    aucs <- apply(query_res, 1, function(y) MESS::auc(query_ns, y))
    aucs <- sort(aucs, decreasing=TRUE)
    query_res <- query_res[names(aucs), ]

    return(query_res)
}

#---------------------


# Bind list of data.frames
#
# rbindlist is fast but data.tables have a different syntax/structure than
# data.frames. Function thus uses rbindlist and returns data.frame with
# restored rownames.
#
# @param l list of data.frames
#
# @return data.frame composed of vertically stacked data.frames in l.

clean_rbindlist <- function(l) {

    #get row names
    l_rnames <- unlist(lapply(l, row.names))

    #bind dfs
    df <- data.table::rbindlist(l)

    #clean up
    class(df) <- "data.frame"
    row.names(df) <- l_rnames

    return(df)
}
