library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

data("sig2Fun_result")

test_that("heatPlot basic output structure", {
    result <- heatPlot(sig2Fun_result)
    expect_type(result, "list")
    expect_named(result, c("heatPlot", "tableHeatPlot"))
    expect_s3_class(result$heatPlot, "ggplot")
    expect_s3_class(result$tableHeatPlot, "data.frame")
})

test_that("heatPlot works with numeric showCategory", {
    result_num <- heatPlot(sig2Fun_result, showCategory = 5)
    expect_s3_class(result_num$heatPlot, "ggplot")
    expect_lte(nrow(result_num$tableHeatPlot), 5 * length(unique(result_num$tableHeatPlot$Gene)))
})

test_that("heatPlot works with specific pathways", {
    gsea <- .extractDF(sig2Fun_result, type = "gseaReadable")
    pathways <- head(gsea@result$ID, 3)
    result_specific <- heatPlot(sig2Fun_result, showCategory = pathways)
    expect_s3_class(result_specific$heatPlot, "ggplot")
    expect_true(all(result_specific$tableHeatPlot$original_name %in% pathways))
})

test_that("heatPlot handles size parameter correctly", {
    res_ok <- heatPlot(sig2Fun_result, size = "pvalue")
    expect_s3_class(res_ok$heatPlot, "ggplot")

    expect_error(
        heatPlot(sig2Fun_result, size = "p.adjust"),
        regexp = "subset|exist"
    )
})

test_that("heatPlot handles color parameter correctly", {
    res_ok <- heatPlot(sig2Fun_result, color = "cor")
    expect_s3_class(res_ok$heatPlot, "ggplot")

    expect_error(
        heatPlot(sig2Fun_result, color = "log2FC"),
        regexp = "subset|exist"
    )
})

test_that("heatPlot preserves factor gene order", {
    result <- heatPlot(sig2Fun_result)
    gene_levels <- levels(result$heatPlot$data$Gene)
    expect_true(length(gene_levels) > 1)
})

test_that("heatPlot produces valid neg_log_size values", {
    result <- heatPlot(sig2Fun_result)
    tbl <- result$tableHeatPlot
    expect_true(all(tbl$neg_log_size >= 0))
})

test_that("heatPlot returns plot and table without error", {
    result <- heatPlot(sig2Fun_result)
    tmp <- tempfile(fileext = ".png")
    suppressWarnings(
        suppressMessages(
            ggsave(tmp, result$heatPlot, width = 6, height = 4)
        )
    )
    expect_true(file.exists(tmp))
    file.remove(tmp)
})
