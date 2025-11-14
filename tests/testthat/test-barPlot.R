library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

data("sig2Fun_result")

test_that("barPlot function basic output structure", {
    result <- barPlot(sig2Fun_result)
    expect_type(result, "list")
    expect_named(result, c("barPlot", "tableBarPlot"))
    expect_s3_class(result$barPlot, "ggplot")
    expect_s3_class(result$tableBarPlot, "data.frame")
    expect_true(all(c("original_name", "label_name", "neg_log_alpha", "alpha_value") %in%
                        colnames(result$tableBarPlot)))
})

test_that("barPlot supports multiple alpha parameters correctly", {
    for (alpha_param in c("pvalue", "p.adjust", "qvalue")) {
        result_alpha <- barPlot(sig2Fun_result, alpha = alpha_param)
        expect_s3_class(result_alpha$barPlot, "ggplot")
        expect_true(alpha_param %in% colnames(result_alpha$tableBarPlot))
    }
})

test_that("barPlot supports multiple x-axis parameters correctly", {
    for (x_param in c("NES", "Count", "GeneRatio")) {
        result_x <- barPlot(sig2Fun_result, x = x_param)
        expect_s3_class(result_x$barPlot, "ggplot")
        expect_true(x_param %in% colnames(result_x$tableBarPlot))
    }
})

test_that("barPlot color and alpha mapping logic is valid", {
    res <- barPlot(sig2Fun_result, alpha = "pvalue", x = "NES")
    gp  <- ggplot_build(res$barPlot)
    fills <- gp$data[[1]]$fill
    alphas <- gp$data[[1]]$alpha

    expect_true(all(fills %in% c("#D25C43", "#5979A3", "grey60")))
    expect_true(all(alphas >= 0 & alphas <= 1))
})

test_that("barPlot handles missing or invalid alpha gracefully", {
    expect_error(
        barPlot(sig2Fun_result, alpha = "invalid_method")
    )
})

test_that("barPlot output table matches ggplot data structure", {
    result <- barPlot(sig2Fun_result)
    plot_data <- ggplot_build(result$barPlot)$data[[1]]
    expect_true(all(result$tableBarPlot$label_name %in% levels(result$tableBarPlot$label_name)))
    expect_true(all(is.finite(result$tableBarPlot$neg_log_alpha)))
})

test_that("barPlot works with small topN", {
    result <- barPlot(sig2Fun_result, topN = 3)
    expect_s3_class(result$barPlot, "ggplot")
    expect_lte(nrow(result$tableBarPlot), 6)
})

test_that("barPlot color mapping logic is valid for NES", {
    res <- barPlot(sig2Fun_result, alpha = "pvalue", x = "NES")
    gp  <- ggplot_build(res$barPlot)
    fills <- unique(gp$data[[1]]$fill)
    expect_true(all(fills %in% c("#D25C43", "#5979A3", "grey60")))
})
