library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

data("demo_GSE181574")

test_that("treePlot function works correctly", {

    result <- treePlot(GSE181574.sigfun)
    expect_type(result, "list")
    expect_named(result, c("treePlot", "Table_treePlot"))
    expect_s3_class(result$treePlot, "ggplot")
    expect_s3_class(result$Table_treePlot, "data.frame")

    result_numeric <- treePlot(GSE181574.sigfun, showCategory=5)
    expect_s3_class(result_numeric$treePlot, "ggplot")
    expect_lte(nrow(result_numeric$Table_treePlot), 5)

    specific_pathway <- c("REACTOME_CELL_CYCLE_MITOTIC",
                          "GOBP_CHROMOSOME_SEPARATION", "REACTOME_SYNTHESIS_OF_DNA",
                          "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")
    result_specific <- treePlot(GSE181574.sigfun, showCategory=specific_pathway)
    expect_s3_class(result_specific$treePlot, "ggplot")
    expect_lte(nrow(result_numeric$Table_treePlot), length(specific_pathway))

})

test_that("treePlot handles errors correctly", {

    expect_snapshot(
        error = TRUE,
        treePlot(GSE181574.sigfun, showCategory=2)
    )

})

test_that("treePlot handles errors correctly", {

    expect_snapshot(
        error = TRUE,
        treePlot(GSE181574.sigfun, showCategory=c(
            "GOMF_STEROID_HYDROXYLASE_ACTIVITY",
            "WP_FATTY_ACID_OMEGAOXIDATION"))
    )

})
