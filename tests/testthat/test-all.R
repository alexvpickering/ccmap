library(ccmap)
library(testthat)



# Toy Data --------------

# create fake drug signatures
genes <- paste("GENE", 1:1000, sep = "_")

set.seed(0)
drug_info <- data.frame(row.names = genes,
                        drug1  = rnorm(1000, sd = 2),
                        drug2  = rnorm(1000, sd = 2),
                        drug3  = rnorm(1000, sd = 2))

# query signature is drug3
query_sig <- drug_info$drug3
names(query_sig) <- genes



# Test Drug Queries -------------------

test_that("Drug queries identify drug3 as perfect match", {

    # best match
    top_drugs <- query_drugs(query_sig, as.matrix(drug_info))
    expect_equal(names(top_drugs)[1], "drug3")

    # match is perfect
    expect_equal(unname(top_drugs[1]), 1)
})

test_that("Drug queries identify reverse drug3 as perfect negative match", {

    # best match
    top_drugs <- query_drugs(-query_sig, as.matrix(drug_info))
    expect_equal(names(top_drugs)[3], "drug3")

    # match is perfect
    expect_equal(unname(top_drugs[3]), -1)
})

test_that("Matches at higher ranks contribute more to overlap", {

    # reverse order of all down-regulated genes
    query_sig  <- sort(query_sig)
    query_sigr <- query_sig
    dn <- query_sigr < 0
    names(query_sigr)[dn] <- rev(names(query_sigr)[dn])

    # signs still match
    expect_equal(sign(query_sig[names(query_sigr)]), sign(query_sigr))

    # overlap is less
    top_drugsr <- query_drugs(query_sigr, as.matrix(drug_info))
    top_drugs  <- query_drugs(query_sig, as.matrix(drug_info))
    expect_lt(top_drugsr[1], top_drugs[1])

    # also reverse order of all up-regulated genes
    up <- query_sigr >= 0
    names(query_sigr)[up] <- rev(names(query_sigr)[up])

    # signs still match
    expect_equal(sign(query_sig[names(query_sigr)]), sign(query_sigr))

    # overlap is even less
    top_drugsrr <- query_drugs(query_sigr, as.matrix(drug_info))
    expect_lt(top_drugsrr[1], top_drugsr[1])
})


# Test Reshaping of Data -------------------

# cmap_es data
cmap_es <- data.frame(g1 = c('d1_g1', 'd2_g1', 'd3_g1'),
                      g2 = c('d1_g2', 'd2_g2', 'd3_g2'))

row.names(cmap_es) <- c('d1', 'd2', 'd3')
cmap_es <- as.matrix(cmap_es)

# use d1 for drug1
drug <- 'd1'
other_drugs <- setdiff(row.names(cmap_es), drug)
drug_es <- ccmap:::rep.row(cmap_es[drug, ], length(other_drugs))

# bind together
Xnet <- cbind(drug_es, cmap_es[other_drugs,, drop=FALSE])

test_that("Xnet is setup correctly", {

    # left side should all be from drug 1
    object <- gsub("_g\\d", "", as.vector(Xnet[,1:2]))
    expect_equal(unique(object), 'd1')

    # g1 in columns 1 and 3
    object <- gsub("d\\d_", "", as.vector(Xnet[,c(1, 3)]))
    expect_equal(unique(object), 'g1')

    # rep.row duplicated rows
    object <- unique(drug_es)
    expect_equal(unique(drug_es), drug_es[1,, drop=FALSE])
})

# NNet predictions
preds <- data.frame(g1 = c('cbo1_g1', 'cbo2_g1'),
                    g2 = c('cbo1_g2', 'cbo2_g2'))

preds <- as.matrix(preds)

# test data for xgboost model
Xgb <- matrix(nrow = length(preds), ncol = 5)
colnames(Xgb) <- c("net_preds", "drug1_dprime", "drug2_dprime",
                    "drug1_vardprime", "drug2_vardprime")

Xgb[, 1] <- as.vector(t(preds))

Xnet1 <- t(Xnet[,1:2])
Xnet2 <- t(Xnet[,3:4])
Xgb[, 2:3] <- cbind(as.vector(Xnet1), as.vector(Xnet2))

# add vardprimes
cmap_var <- data.frame(d1 = c('d1_g1','d1_g2'),
                       d2 = c('d2_g1','d2_g2'),
                       d3 = c('d3_g1','d3_g2'))

cmap_var <- as.matrix(cmap_var)

Xgb[, 4] <- cmap_var[, drug]
Xgb[, 5] <- as.vector(cmap_var[, other_drugs])

test_that("Xgb is setup correctly", {

    # each row is for the same gene
    for (i in 1:nrow(Xgb)) {
        object <- gsub("\\w+_", "", as.vector(Xgb[i, ]))
        expect_length(unique(object), 1)
    }

    # dprimes and vardprimes are for same drug and gene
    expect_equal(Xgb[, 'drug1_dprime'], Xgb[, 'drug1_vardprime'])
    expect_equal(Xgb[, 'drug2_dprime'], Xgb[, 'drug2_vardprime'])
})

# setup combos_es
combos_es <- Xgb[, 'net_preds']
dim(combos_es) <- c(2, length(combos_es)/2)
colnames(combos_es) <- paste(drug, other_drugs, sep = " + ")
row.names(combos_es) <- c('g1', 'g2')

test_that("combos_es is setup correctly", {

    # each column is a combo
    for (i in 1:ncol(combos_es)) {
        object <- gsub("_\\w+", "", as.vector(combos_es[, i]))
        expect_length(unique(object), 1)
    }

    # each row is a gene
    for (i in 1:nrow(combos_es)) {
        object <- gsub("\\w+_", "", as.vector(combos_es[i, ]))
        expect_length(unique(object), 1)
    }
})


# Test NNet prediction function ----------


# setup net1
nsamples <- 5
ngenes   <- 10
nhidden  <- 2

Xnet <- matrix(1, nsamples, ngenes*2)

net1 <- list(W1 = matrix(1, ngenes*2, nhidden),
             W2 = matrix(1, nhidden,  ngenes),
             b1 = rep(0, nhidden),
             b2 = rep(0, ngenes))

test_that("predict.nnet produces correct dimension output", {
    preds <- ccmap:::predict.nnet(net1, Xnet)

    expect_equal(dim(preds), c(nsamples, ngenes))
})

# Test get_bins ------

test_that("get_bins is working properly", {

    # if ask for more bins than n, gives n bins of length 1
    bins <- ccmap:::get_bins(10, 11)
    expect_length(bins, 10)

    bin_lengths <- unlist(lapply(bins, length), use.names = FALSE)
    expect_equal(bin_lengths, rep(1, 10))

    # if ask for 1 bin, gives 1 bin of length n
    bins <- ccmap:::get_bins(10, 1)
    expect_length(bins, 1)
    expect_length(bins[[1]], 10)

    # difference between length of bins is never more than 1
    for (i in 1:100) {
        bins <- ccmap:::get_bins(100, i)
        bin_lengths <- unlist(lapply(bins, length), use.names = FALSE)
        expect_lte(max(bin_lengths) - min(bin_lengths), 1)
    }
})


# cleanup -----
rm(drug_info, genes, query_sig, cmap_es, cmap_var, combos_es, drug_es, preds,
   Xgb, Xnet, Xnet1, Xnet2, drug, other_drugs, nsamples, ngenes, nhidden, Xnet,
   net1)
