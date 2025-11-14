library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

data("sig2Fun_result")

test_that("upsetPlot basic output structure", {
    result <- upsetPlot(sig2Fun_result)
    expect_type(result, "list")
    expect_named(result, c("upsetPlot", "tableUpsetPlot"))
    expect_s3_class(result$upsetPlot, "ggplot")
    expect_s3_class(result$tableUpsetPlot, "data.frame")
})

test_that("upsetPlot type parameter works", {
    result_box <- upsetPlot(sig2Fun_result, type = "box")
    result_bar <- upsetPlot(sig2Fun_result, type = "bar")
    expect_s3_class(result_box$upsetPlot, "ggplot")
    expect_s3_class(result_bar$upsetPlot, "ggplot")
})

test_that("upsetPlot showCategory parameter works", {
    result_numeric <- upsetPlot(sig2Fun_result, showCategory = 5)
    expect_s3_class(result_numeric$upsetPlot, "ggplot")

    specific_pathway <- c(
        "REACTOME_CELL_CYCLE_MITOTIC",
        "GOBP_CHROMOSOME_SEPARATION"
    )
    result_specific <- upsetPlot(sig2Fun_result, showCategory = specific_pathway)
    expect_s3_class(result_specific$upsetPlot, "ggplot")
})

test_that("upsetPlot handles invalid parameters", {
    expect_error(
        upsetPlot(sig2Fun_result, type = "invalid_method"),
        regexp = "'arg' should be one of"
    )
})


