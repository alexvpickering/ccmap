#---------------------


#' Title
#'
#' @param scores
#' @param drug_name
#' @param drug_info
#' @param query_ns
#' @param es
#'
#' @return
#' @export
#'
#' @examples
test_ma <- function(scores, drug_name, drug_info, query_ns=2^(4:13), es=TRUE) {

    #get top up/dn drug genes
    drug_genes <- get_drug_genes(drug_info, scores$meta, es=es)

    #final results
    res_cons <- list()
    res_meta <- c()

    for (query_n in query_ns) {

        #perform query of size query_n with meta result
        meta_top <- get_top_drugs(scores$meta, drug_info,
                                  query_n, drug_genes, es=es)

        meta_rank <- get_drug_rank(drug_name, meta_top)
        res_meta <- c(res_meta, meta_rank)

        cons_ranks <- c()
        for (con in scores$contrasts) {

            #perform query of size query_n with each contrast
            con_top <- get_top_drugs(con, drug_info,
                                     query_n, drug_genes, es=es)

            con_rank <- get_drug_rank(drug_name, con_top)

            #combine rank with those of other contrasts
            cons_ranks <- c(cons_ranks, con_rank)
        }
        #store ranks of each contrast for size query_n
        res_cons[[as.character(query_n)]] <- cons_ranks
    }
    #add names to res_meta
    names(res_meta) <- as.character(query_ns)

    ma_res <- list(meta=res_meta, cons=res_cons)
    return(ma_res)
}


#---------------------


#' Title
#'
#' @import ggplot2
#'
#' @param ma_res
#'
#' @return
#' @export
#'
#' @examples
plot_ma_res <- function (ma_res, ymax=1310, ytick=100) {

    cons_res <- data.frame(ma_res$cons)
    meta_res <- t(data.frame(ma_res$meta))
    colnames(cons_res) <- colnames(meta_res)

    cons_res <- reshape2::melt(cons_res)
    meta_res <- reshape2::melt(meta_res)[,-1]
    colnames(cons_res) <- c("query_n", "rank")
    colnames(meta_res) <- c("query_n", "rank")

    meta_res$meta <- 1
    cons_res$meta <- 0

    full_res <- rbind(meta_res, cons_res)
    full_res$rank <- as.integer(full_res$rank)
    full_res$query_n <- as.integer(full_res$query_n)

    p <- ggplot(full_res, aes(query_n, rank))

    p + stat_summary(fun.y="median", colour = "red",
                     size=10, geom="point", shape=45) +

        geom_jitter(width=0.3, aes(alpha=factor(meta))) +

        scale_y_continuous(breaks=seq(0, 1310, ytick),
                           limits=c(0, 1310)) +

        scale_x_continuous(breaks=2^(4:13),
                           limits=c(12, 10000),
                           trans=scales::log2_trans()) +

        coord_cartesian(ylim=c(0, ymax))
}
