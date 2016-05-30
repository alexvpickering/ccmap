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

test_that("Drug queries identify drug3 as best match", {

    top_drugs <- query_drugs(query_sig, as.matrix(drug_info))
    expect_equal(row.names(top_drugs)[1], "drug3")
})

# cleanup
rm(drug_info, genes, query_sig)
