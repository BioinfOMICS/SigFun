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

test_that("lollipopPlot numeric showCategory with NES", {
    # When x = "NES", showCategory = N gives up to 2N pathways (N positive + N negative)
    result_num <- lollipopPlot(sig2Fun_result, showCategory = 5, x = "NES")
    expect_s3_class(result_num$lollipopPlot, "ggplot")
    # Should get up to 10 pathways (5 positive + 5 negative)
    expect_lte(nrow(result_num$tableLollipopPlot), 10)

    # Check that we have both positive and negative NES if they exist
    if (nrow(result_num$tableLollipopPlot) > 1) {
        nes_values <- result_num$tableLollipopPlot$NES
        # If both positive and negative exist, check we have both
        gseaRaw <- .extractDF(sig2Fun_result, type = "gseaRaw")
        has_positive <- any(gseaRaw@result$NES > 0)
        has_negative <- any(gseaRaw@result$NES < 0)
        if (has_positive && has_negative) {
            expect_true(any(nes_values > 0) && any(nes_values < 0))
        }
    }
})

test_that("lollipopPlot numeric showCategory with non-NES", {
    # When x != "NES", showCategory = N gives exactly N pathways
    result_count <- lollipopPlot(sig2Fun_result, showCategory = 5, x = "Count")
    expect_s3_class(result_count$lollipopPlot, "ggplot")
    expect_lte(nrow(result_count$tableLollipopPlot), 5)

    result_ratio <- lollipopPlot(sig2Fun_result, showCategory = 5, x = "GeneRatio")
    expect_s3_class(result_ratio$lollipopPlot, "ggplot")
    expect_lte(nrow(result_ratio$tableLollipopPlot), 5)
})

test_that("lollipopPlot specific showCategory", {
    gseaRaw <- .extractDF(sig2Fun_result, type = "gseaRaw")
    pathways <- head(gseaRaw@result$ID, 3)
    result_specific <- lollipopPlot(sig2Fun_result, showCategory = pathways)
    expect_s3_class(result_specific$lollipopPlot, "ggplot")
    expect_equal(nrow(result_specific$tableLollipopPlot), 3)
    expect_true(all(result_specific$tableLollipopPlot$original_name %in% pathways))
})

test_that("lollipopPlot handles invalid parameters", {
    expect_error(lollipopPlot(sig2Fun_result, x = "invalid"), regexp = "should be one of")
    expect_error(lollipopPlot(sig2Fun_result, alpha = "invalid"), regexp = "should be one of")
    expect_error(lollipopPlot(sig2Fun_result, size = "invalid"), regexp = "should be one of")
})

test_that("lollipopPlot color mapping for NES", {
    res <- lollipopPlot(sig2Fun_result, x = "NES", alpha = "pvalue")
    gp <- ggplot_build(res$lollipopPlot)

    # When x = "NES", colors should be red (#b3200a) and blue (#08306b)
    # Find the geom_point layer (should be the last one if addSeg = TRUE)
    point_layer <- tail(gp$data, 1)[[1]]
    expect_true(all(point_layer$colour %in% c("#b3200a", "#08306b")))
    expect_true(all(point_layer$alpha >= 0 & point_layer$alpha <= 1))
})

test_that("lollipopPlot color mapping for non-NES", {
    res <- lollipopPlot(sig2Fun_result, x = "Count", alpha = "pvalue")
    gp <- ggplot_build(res$lollipopPlot)

    # When x != "NES", all colors should be red (#b3200a)
    point_layer <- tail(gp$data, 1)[[1]]
    expect_true(all(point_layer$colour == "#b3200a"))
})

test_that("lollipopPlot addSeg parameter works", {
    res_seg <- lollipopPlot(sig2Fun_result, addSeg = TRUE)
    res_noseg <- lollipopPlot(sig2Fun_result, addSeg = FALSE)

    expect_s3_class(res_seg$lollipopPlot, "ggplot")
    expect_s3_class(res_noseg$lollipopPlot, "ggplot")

    # Check number of layers
    gp_seg <- ggplot_build(res_seg$lollipopPlot)
    gp_noseg <- ggplot_build(res_noseg$lollipopPlot)

    # With segment should have one more layer than without
    expect_equal(length(gp_seg$data), length(gp_noseg$data) + 1)
})

test_that("lollipopPlot lineType and lineSize parameters work", {
    res1 <- lollipopPlot(sig2Fun_result, lineType = "dashed", lineSize = 0.5, addSeg = TRUE)
    res2 <- lollipopPlot(sig2Fun_result, lineType = "solid", lineSize = 2, addSeg = TRUE)

    expect_s3_class(res1$lollipopPlot, "ggplot")
    expect_s3_class(res2$lollipopPlot, "ggplot")

    # Check that segments exist and have correct properties
    gp1 <- ggplot_build(res1$lollipopPlot)
    gp2 <- ggplot_build(res2$lollipopPlot)

    # First layer should be segments
    expect_true("GeomSegment" %in% class(res1$lollipopPlot$layers[[1]]$geom))
    expect_true("GeomSegment" %in% class(res2$lollipopPlot$layers[[1]]$geom))
})

test_that("lollipopPlot output table has correct columns", {
    result <- lollipopPlot(sig2Fun_result)
    expect_true(all(c("original_name", "Description", "neg_log_alpha") %in%
                        colnames(result$tableLollipopPlot)))

    # Check that dynamic columns exist (based on x, alpha, size parameters)
    expect_true("NES" %in% colnames(result$tableLollipopPlot))  # default x
    expect_true("p.adjust" %in% colnames(result$tableLollipopPlot))  # default alpha
    expect_true("Count" %in% colnames(result$tableLollipopPlot))  # default size
})

test_that("lollipopPlot dynamic columns in output table", {
    result1 <- lollipopPlot(sig2Fun_result, x = "GeneRatio", alpha = "pvalue", size = "GeneRatio")
    expect_true(all(c("GeneRatio", "pvalue") %in% colnames(result1$tableLollipopPlot)))

    result2 <- lollipopPlot(sig2Fun_result, x = "Count", alpha = "qvalue", size = "Count")
    expect_true(all(c("Count", "qvalue") %in% colnames(result2$tableLollipopPlot)))
})

test_that("lollipopPlot handles fontSize parameter", {
    res_small <- lollipopPlot(sig2Fun_result, fontSize = 8)
    res_large <- lollipopPlot(sig2Fun_result, fontSize = 14)

    expect_s3_class(res_small$lollipopPlot, "ggplot")
    expect_s3_class(res_large$lollipopPlot, "ggplot")
})

test_that("lollipopPlot works with all x-axis options", {
    for (x_param in c("NES", "Count", "GeneRatio")) {
        result <- lollipopPlot(sig2Fun_result, x = x_param)
        expect_s3_class(result$lollipopPlot, "ggplot")
        expect_true(x_param %in% colnames(result$tableLollipopPlot))
    }
})

test_that("lollipopPlot breaklineN parameter works", {
    res_short <- lollipopPlot(sig2Fun_result, breaklineN = 20)
    res_long <- lollipopPlot(sig2Fun_result, breaklineN = 50)

    expect_s3_class(res_short$lollipopPlot, "ggplot")
    expect_s3_class(res_long$lollipopPlot, "ggplot")
})

test_that("lollipopPlot saves plot without error", {
    result <- lollipopPlot(sig2Fun_result)
    tmp <- tempfile(fileext = ".png")
    suppressWarnings(
        suppressMessages(
            ggsave(tmp, result$lollipopPlot, width = 6, height = 4)
        )
    )
    expect_true(file.exists(tmp))
    file.remove(tmp)
})

test_that("lollipopPlot x-axis range is correct", {
    # For NES, x-axis should include 0 and have buffer
    res_nes <- lollipopPlot(sig2Fun_result, x = "NES")
    gp_nes <- ggplot_build(res_nes$lollipopPlot)
    x_range_nes <- gp_nes$layout$panel_scales_x[[1]]$range$range
    expect_true(x_range_nes[1] <= 0 && x_range_nes[2] >= 0)

    # For Count/GeneRatio, x-axis should start from 0
    res_count <- lollipopPlot(sig2Fun_result, x = "Count")
    gp_count <- ggplot_build(res_count$lollipopPlot)
    x_range_count <- gp_count$layout$panel_scales_x[[1]]$range$range
    expect_equal(x_range_count[1], 0)
})
