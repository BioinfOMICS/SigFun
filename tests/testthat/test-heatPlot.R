library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)
data("sig2Fun_result")

test_that("lollipopPlot basic output structure", {
    result <- lollipopPlot(sig2Fun_result)
    expect_type(result, "list")
    expect_named(result, c("lollipopPlot", "tableLollipopPlot"))
    expect_s3_class(result$lollipopPlot, "ggplot")
    expect_s3_class(result$tableLollipopPlot, "data.frame")
})

test_that("lollipopPlot numeric and specific showCategory", {
    result_num <- lollipopPlot(sig2Fun_result, showCategory = 5)
    expect_s3_class(result_num$lollipopPlot, "ggplot")
    expect_lte(nrow(result_num$tableLollipopPlot), 5)
    gseaRaw <- .extractDF(sig2Fun_result, type = "gseaRaw")
    pathways <- head(gseaRaw@result$ID, 3)
    result_specific <- lollipopPlot(sig2Fun_result, showCategory = pathways)
    expect_s3_class(result_specific$lollipopPlot, "ggplot")
    expect_lte(nrow(result_specific$tableLollipopPlot), 3)
})

test_that("lollipopPlot handles invalid parameters", {
    expect_error(lollipopPlot(sig2Fun_result, x = "invalid"), regexp = "should be one of")
    expect_error(lollipopPlot(sig2Fun_result, alpha = "invalid"), regexp = "should be one of")
    expect_error(lollipopPlot(sig2Fun_result, size = "invalid"), regexp = "should be one of")
})

test_that("heatPlot color scale and fill values behave as expected", {
    res <- heatPlot(sig2Fun_result, color = "cor", size = "pvalue")
    gp <- ggplot_build(res$heatPlot)
    fill_colors <- gp$data[[1]]$fill
    is_valid_color <- grepl("^#[0-9A-Fa-f]{6}$", fill_colors) | fill_colors == "white"
    expect_true(all(is_valid_color))
    expect_gt(length(unique(fill_colors)), 1)
    expect_true(all(gp$data[[1]]$size >= 0))
    expect_true(all(is.finite(gp$data[[1]]$size)))
})

test_that("lollipopPlot addSeg parameter works", {
    res_seg <- lollipopPlot(sig2Fun_result, addSeg = TRUE)
    res_noseg <- lollipopPlot(sig2Fun_result, addSeg = FALSE)
    expect_s3_class(res_seg$lollipopPlot, "ggplot")
    expect_s3_class(res_noseg$lollipopPlot, "ggplot")
})

test_that("lollipopPlot lineType and lineSize parameters work", {
    res1 <- lollipopPlot(sig2Fun_result, lineType = "dashed", lineSize = 0.5)
    res2 <- lollipopPlot(sig2Fun_result, lineType = "solid", lineSize = 2)
    expect_s3_class(res1$lollipopPlot, "ggplot")
    expect_s3_class(res2$lollipopPlot, "ggplot")
})

test_that("lollipopPlot output table has correct columns", {
    result <- lollipopPlot(sig2Fun_result)
    expect_true(all(c("original_name", "Description", "neg_log_alpha") %in%
                        colnames(result$tableLollipopPlot)))
})

test_that("lollipopPlot handles fontSize parameter", {
    res_small <- lollipopPlot(sig2Fun_result, fontSize = 8)
    res_large <- lollipopPlot(sig2Fun_result, fontSize = 14)
    expect_s3_class(res_small$lollipopPlot, "ggplot")
    expect_s3_class(res_large$lollipopPlot, "ggplot")
})

test_that("lollipopPlot works with non-NES x-axis", {
    for (x_param in c("Count", "GeneRatio")) {
        result <- lollipopPlot(sig2Fun_result, x = x_param)
        expect_s3_class(result$lollipopPlot, "ggplot")
    }
})

test_that("lollipopPlot saves plot without error", {
    result <- lollipopPlot(sig2Fun_result)
    temp_file <- tempfile(fileext = ".png")
    suppressWarnings(
        suppressMessages(
            ggsave(temp_file, result$lollipopPlot, width = 6, height = 4)
            )
        )
    expect_true(file.exists(temp_file))
    expect_gt(file.size(temp_file), 0)
    file.remove(temp_file)
})


