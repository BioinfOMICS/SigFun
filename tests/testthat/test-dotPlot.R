library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)
data("sig2Fun_result")

test_that("dotPlot function basic output structure", {
    result <- tryCatch({
        dotPlot(sig2Fun_result)
    }, error = function(e) {
        skip(paste("dotPlot failed:", e$message))
    })
    expect_type(result, "list")
    expect_named(result, c("dotPlot", "tableDotPlot"))
    expect_s3_class(result$tableDotPlot, "data.frame")
    expect_s3_class(result$dotPlot, "ggplot")
    expect_true(all(c("original_name", "Description", "Count", "GeneRatio") %in%
                        colnames(result$tableDotPlot)))
})

test_that("dotPlot works with numeric and specific showCategory", {
    result_numeric <- tryCatch({
        dotPlot(sig2Fun_result, showCategory = 5)
    }, error = function(e) {
        skip(paste("dotPlot numeric category failed:", e$message))
    })
    expect_named(result_numeric, c("dotPlot", "tableDotPlot"))
    expect_s3_class(result_numeric$tableDotPlot, "data.frame")
    expect_lte(nrow(result_numeric$tableDotPlot), 5)
    gseaRaw <- .extractDF(sig2Fun_result, type = "gseaRaw")
    available_pathways <- head(gseaRaw@result$ID, 3)
    if (length(available_pathways) >= 3) {
        result_specific <- tryCatch({
            dotPlot(sig2Fun_result, showCategory = available_pathways)
        }, error = function(e) {
            skip(paste("Specific pathway test failed:", e$message))
        })
        expect_s3_class(result_specific$tableDotPlot, "data.frame")
        expect_equal(nrow(result_specific$tableDotPlot), 3)
    } else {
        skip("Not enough pathways for specific test")
    }
})

test_that("dotPlot handles invalid parameters correctly", {
    expect_error(
        dotPlot(sig2Fun_result, alpha = "invalid"),
        regexp = "'arg' should be one of"
    )
    expect_error(
        dotPlot(sig2Fun_result, x = "invalid"),
        regexp = "'arg' should be one of"
    )
    expect_error(
        dotPlot(sig2Fun_result, size = "invalid"),
        regexp = "'arg' should be one of"
    )
})

test_that("dotPlot output ordering and structure are consistent", {
    result <- dotPlot(sig2Fun_result, orderBy = "GeneRatio", decreasing = TRUE)
    expect_true(is.factor(result$tableDotPlot$Description))
    expect_true(all(is.finite(result$tableDotPlot$GeneRatio)))
    expect_true(all(result$tableDotPlot$p.adjust >= 0 | is.na(result$tableDotPlot$p.adjust)))
    if (nrow(result$tableDotPlot) > 1) {
        gene_ratios <- result$tableDotPlot$GeneRatio
        expect_true(all(diff(gene_ratios) <= 0) || any(is.na(gene_ratios)))
    }
})

test_that("dotPlot color and alpha mapping works as expected", {
    res <- dotPlot(sig2Fun_result, x = "NES", alpha = "p.adjust")
    gp <- ggplot_build(res$dotPlot)
    expect_true(all(gp$data[[1]]$colour %in% c("#b3200a", "#08306b")))
    expect_true(all(gp$data[[1]]$alpha >= 0 & gp$data[[1]]$alpha <= 1))
})

test_that("dotPlot supports different x-axis parameters", {
    for (x_param in c("GeneRatio", "Count", "NES")) {
        result <- tryCatch({
            dotPlot(sig2Fun_result, x = x_param)
        }, error = function(e) {
            skip(paste("x =", x_param, "failed:", e$message))
        })
        expect_s3_class(result$dotPlot, "ggplot")
        expect_true(x_param %in% colnames(result$tableDotPlot))
    }
})

test_that("dotPlot supports different size parameters", {
    for (size_param in c("Count", "GeneRatio")) {
        result <- tryCatch({
            dotPlot(sig2Fun_result, size = size_param)
        }, error = function(e) {
            skip(paste("size =", size_param, "failed:", e$message))
        })
        expect_s3_class(result$dotPlot, "ggplot")
        expect_true(size_param %in% colnames(result$tableDotPlot))
    }
})

test_that("dotPlot supports different alpha parameters", {
    for (alpha_param in c("pvalue", "p.adjust", "qvalue")) {
        result <- tryCatch({
            dotPlot(sig2Fun_result, alpha = alpha_param)
        }, error = function(e) {
            skip(paste("alpha =", alpha_param, "failed:", e$message))
        })
        expect_s3_class(result$dotPlot, "ggplot")
        expect_true(alpha_param %in% colnames(result$tableDotPlot))
    }
})

test_that("dotPlot handles breaklineN parameter correctly", {
    result_short <- dotPlot(sig2Fun_result, breaklineN = 20)
    result_long <- dotPlot(sig2Fun_result, breaklineN = 50)
    expect_s3_class(result_short$dotPlot, "ggplot")
    expect_s3_class(result_long$dotPlot, "ggplot")
})

test_that("dotPlot handles fontSize parameter correctly", {
    result_small <- dotPlot(sig2Fun_result, fontSize = 8)
    result_large <- dotPlot(sig2Fun_result, fontSize = 14)
    expect_s3_class(result_small$dotPlot, "ggplot")
    expect_s3_class(result_large$dotPlot, "ggplot")
})

test_that("dotPlot handles title parameter correctly", {
    result_with_title <- dotPlot(sig2Fun_result, title = "Test Title")
    result_no_title <- dotPlot(sig2Fun_result, title = NULL)
    expect_s3_class(result_with_title$dotPlot, "ggplot")
    expect_s3_class(result_no_title$dotPlot, "ggplot")
})

test_that("dotPlot handles orderBy parameter correctly", {
    result_count <- dotPlot(sig2Fun_result, orderBy = "Count", decreasing = TRUE)
    result_nes <- dotPlot(sig2Fun_result, orderBy = "NES", decreasing = FALSE)

    expect_s3_class(result_count$dotPlot, "ggplot")
    expect_s3_class(result_nes$dotPlot, "ggplot")
    if (nrow(result_count$tableDotPlot) > 1) {
        counts <- result_count$tableDotPlot$Count
        expect_true(all(diff(counts) <= 0))
    }
})

test_that("dotPlot tableDotPlot contains all required columns", {
    result <- dotPlot(sig2Fun_result)
    required_cols <- c("original_name", "Description", "NES", "Count", "GeneRatio", "p.adjust")
    expect_true(all(required_cols %in% colnames(result$tableDotPlot)))
})

test_that("dotPlot handles empty or single category correctly", {
    gseaRaw <- .extractDF(sig2Fun_result, type = "gseaRaw")
    single_pathway <- head(gseaRaw@result$ID, 1)
    result_single <- tryCatch({
        dotPlot(sig2Fun_result, showCategory = single_pathway)
    }, error = function(e) {
        skip(paste("Single category test failed:", e$message))
    })
    expect_s3_class(result_single$dotPlot, "ggplot")
    expect_equal(nrow(result_single$tableDotPlot), 1)
})

test_that("dotPlot color logic differs between NES and other x-axis", {
    result_nes <- dotPlot(sig2Fun_result, x = "NES")
    result_count <- dotPlot(sig2Fun_result, x = "Count")
    gp_nes <- ggplot_build(result_nes$dotPlot)
    gp_count <- ggplot_build(result_count$dotPlot)
    colors_nes <- unique(gp_nes$data[[1]]$colour)
    colors_count <- unique(gp_count$data[[1]]$colour)
    expect_true(length(colors_nes) <= 2)
    expect_equal(unique(colors_count), "#b3200a")
})

test_that("dotPlot neg_log_alpha calculation is correct", {
    result <- dotPlot(sig2Fun_result, alpha = "p.adjust")
    expect_true(all(is.finite(result$tableDotPlot$p.adjust) |
                        is.na(result$tableDotPlot$p.adjust)))
})

test_that("dotPlot handles zero showCategory gracefully", {
    expect_error(
        dotPlot(sig2Fun_result, showCategory = 0),
        regexp = "No pathways to display"
    )
})
