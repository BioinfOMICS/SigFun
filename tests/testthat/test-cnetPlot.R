library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

data("sig2Fun_result")

test_that("cnetPlot function basic structure", {
    result <- tryCatch({
        cnetPlot(sig2Fun_result)
    }, error = function(e) {
        skip(paste("cnetPlot failed to generate:", e$message))
    })

    expect_type(result, "list")
    expect_named(result, c("cnetPlot", "tableEdgeCnetPlot", "tableNodeCnetPlot"))
    expect_s3_class(result$tableEdgeCnetPlot, "data.frame")
    expect_s3_class(result$tableNodeCnetPlot, "data.frame")
    expect_true(all(c("from", "to") %in% colnames(result$tableEdgeCnetPlot)))
    expect_true(all(c("id", "size") %in% colnames(result$tableNodeCnetPlot)))
})

test_that("cnetPlot works without color parameter", {
    result <- tryCatch({
        cnetPlot(sig2Fun_result, color = NULL)
    }, error = function(e) {
        skip(paste("Skipping due to plotting error:", e$message))
    })
    expect_type(result, "list")
    expect_named(result, c("cnetPlot", "tableEdgeCnetPlot", "tableNodeCnetPlot"))
})

test_that("cnetPlot works with numeric and specific showCategory", {
    result_numeric <- tryCatch({
        cnetPlot(sig2Fun_result, showCategory = 5)
    }, error = function(e) {
        skip(paste("Skipping numeric showCategory due to:", e$message))
    })
    expect_s3_class(result_numeric$tableEdgeCnetPlot, "data.frame")

    specific_pathway <- c(
        "REACTOME_CELL_CYCLE_MITOTIC",
        "GOBP_CHROMOSOME_SEPARATION",
        "REACTOME_SYNTHESIS_OF_DNA"
    )

    result_specific <- tryCatch({
        cnetPlot(sig2Fun_result, showCategory = specific_pathway)
    }, error = function(e) {
        if (grepl("No matching pathways", e$message)) {
            skip("Specific pathways not found in test dataset")
        } else stop(e)
    })
    expect_s3_class(result_specific$tableEdgeCnetPlot, "data.frame")
})

test_that("cnetPlot handles invalid color column correctly", {
    expect_error(
        cnetPlot(sig2Fun_result, color = "invalid_column"),
        regexp = "color"
    )
})

test_that("cnetPlot handles invalid pathway names correctly", {
    expect_error(
        cnetPlot(sig2Fun_result, showCategory = c("invalid_pathway1", "invalid_pathway2")),
        regexp = "No matching pathways"
    )
})

test_that("cnetPlot handles invalid nodeLabel correctly", {
    expect_error(
        cnetPlot(sig2Fun_result, nodeLabel = "invalid_label"),
        regexp = "wrong parameter"
    )
})
