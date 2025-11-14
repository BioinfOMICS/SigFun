library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)
data("sig2Fun_result")

test_that("ridgePlot basic output structure", {
    result <- ridgePlot(sig2Fun_result)
    expect_type(result, "list")
    expect_named(result, c("ridgePlot", "tableRidgePlot"))
    expect_s3_class(result$ridgePlot, "ggplot")
    expect_s3_class(result$tableRidgePlot, "data.frame")
})

test_that("ridgePlot numeric and specific showCategory", {
    result_num <- ridgePlot(sig2Fun_result, showCategory = 5)
    expect_s3_class(result_num$ridgePlot, "ggplot")
    gseaRaw <- .extractDF(sig2Fun_result, type = "gseaRaw")
    pathways <- head(gseaRaw@result$ID, 3)
    result_specific <- ridgePlot(sig2Fun_result, showCategory = pathways)
    expect_s3_class(result_specific$ridgePlot, "ggplot")
})

test_that("ridgePlot handles invalid parameters", {
    expect_error(ridgePlot(sig2Fun_result, fill = "invalid"), regexp = "should be one of")
})

test_that("ridgePlot output table has correct columns", {
    result <- ridgePlot(sig2Fun_result)
    expect_true(all(c(
        "original_name", "label_name", "raw_fill",
        "neg_log10_fill", "ranking_metric", "median_ranking_metric"
    ) %in% colnames(result$tableRidgePlot)))
})

test_that("ridgePlot decreasing parameter works", {
    res1 <- ridgePlot(sig2Fun_result, decreasing = TRUE)
    res2 <- ridgePlot(sig2Fun_result, decreasing = FALSE)
    expect_s3_class(res1$ridgePlot, "ggplot")
    expect_s3_class(res2$ridgePlot, "ggplot")
})

test_that("ridgePlot breaklineN and fontSize work", {
    res_small <- ridgePlot(sig2Fun_result, breaklineN = 20, fontSize = 8)
    res_large <- ridgePlot(sig2Fun_result, breaklineN = 50, fontSize = 14)
    expect_s3_class(res_small$ridgePlot, "ggplot")
    expect_s3_class(res_large$ridgePlot, "ggplot")
})

test_that("ridgePlot color scale uses fill column correctly", {
    res <- ridgePlot(sig2Fun_result, fill = "p.adjust")
    gp <- ggplot_build(res$ridgePlot)
    expect_true(all(gp$plot$labels$fill %in% c("-log10(p.adjust)")))
})

test_that("ridgePlot saves plot without error", {
    result <- ridgePlot(sig2Fun_result)
    tmp <- tempfile(fileext = ".png")
    suppressWarnings(
        suppressMessages(
            ggsave(tmp, result$ridgePlot, width = 6, height = 4)
        )
    )
    expect_true(file.exists(tmp))
    file.remove(tmp)
})
