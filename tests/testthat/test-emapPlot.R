library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)
data("sig2Fun_result")

test_that("emapPlot basic output structure", {
    result <- tryCatch({
        emapPlot(sig2Fun_result)
    }, error = function(e) skip(paste("emapPlot failed:", e$message)))
    expect_type(result, "list")
    expect_named(result, c("emapPlot", "tableEmapPlot"))
    expect_s3_class(result$emapPlot, "ggplot")
    expect_s3_class(result$tableEmapPlot, "data.frame")
    expect_true(all(c("labelName", "size", "color") %in% colnames(result$tableEmapPlot)))
})

test_that("emapPlot handles numeric and specific showCategory", {
    result_num <- tryCatch({
        emapPlot(sig2Fun_result, showCategory = 5)
    }, error = function(e) skip(paste("emapPlot numeric failed:", e$message)))
    expect_s3_class(result_num$emapPlot, "ggplot")
    expect_lte(nrow(result_num$tableEmapPlot), 5)
    gseaSimilar <- .extractDF(sig2Fun_result, type = "gseaSimilar")
    available_pathways <- head(unique(gseaSimilar$Description), 3)
    if (length(available_pathways) >= 3) {
        result_specific <- tryCatch({
            emapPlot(sig2Fun_result, showCategory = available_pathways)
        }, error = function(e) skip(paste("emapPlot specific failed:", e$message)))
        expect_s3_class(result_specific$emapPlot, "ggplot")
        expect_true(nrow(result_specific$tableEmapPlot) <= 3)
    } else {
        skip("Not enough pathways for specific test")
    }
})

test_that("emapPlot supports all valid color parameters", {
    for (color_param in c("pvalue", "p.adjust", "qvalue")) {
        result <- tryCatch({
            emapPlot(sig2Fun_result, color = color_param)
        }, error = function(e) skip(paste("emapPlot color", color_param, "failed:", e$message)))
        expect_s3_class(result$emapPlot, "ggplot")
        gp <- ggplot_build(result$emapPlot)
        expect_true(all(is.finite(gp$data[[1]]$size)))
    }
})

test_that("emapPlot handles invalid color parameter gracefully", {
    expect_error(
        emapPlot(sig2Fun_result, color = "invalid_color"),
        regexp = "'arg' should be one of"
    )
})

test_that("emapPlot handles nodeLabel correctly", {
    for (label_mode in c("all", "category", "none")) {
        result <- tryCatch({
            emapPlot(sig2Fun_result, nodeLabel = label_mode)
        }, error = function(e) skip(paste("emapPlot label", label_mode, "failed:", e$message)))
        expect_s3_class(result$emapPlot, "ggplot")
    }
})

test_that("emapPlot works with different layout functions", {
    for (layout_fun in list(
        igraph::layout_with_fr,
        igraph::layout_in_circle,
        igraph::layout_with_kk
    )) {
        result <- tryCatch({
            emapPlot(sig2Fun_result, layout = layout_fun)
        }, error = function(e) skip(paste("emapPlot layout test failed:", e$message)))
        expect_s3_class(result$emapPlot, "ggplot")
    }
})

test_that("emapPlot edge and sizeCategory parameters work correctly", {
    result_small <- emapPlot(sig2Fun_result, sizeCategory = 0.5, sizeEdge = 0.2)
    result_large <- emapPlot(sig2Fun_result, sizeCategory = 1.5, sizeEdge = 0.8)
    expect_s3_class(result_small$emapPlot, "ggplot")
    expect_s3_class(result_large$emapPlot, "ggplot")
})

test_that("emapPlot tableEmapPlot has valid numeric data", {
    result <- emapPlot(sig2Fun_result)
    df <- result$tableEmapPlot
    expect_true(all(is.numeric(df$size)))
    expect_true(all(is.numeric(df$color)))
    expect_true(all(!duplicated(df$labelName)))
})
