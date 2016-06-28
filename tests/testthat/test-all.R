library(ccmap)
library(testthat)



# -------------- Toy Data

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



# ------------------- Test Drug Queries

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





# cleanup
rm(drug_info, genes, query_sig)

