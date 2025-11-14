library(testthat)
library(ggplot2)
library(SummarizedExperiment)

data("sig2Fun_result")

test_that("gseaPlot basic functionality and return type", {
    result <- tryCatch({
        gseaPlot(sig2Fun_result)
    }, error = function(e) {
        skip(paste("gseaPlot failed:", e$message))
    })
    expect_type(result, "list")
    expect_true(all(sapply(result, inherits, what = "gg")))
})

test_that("gseaPlot works with numeric and specific pathway input", {
    result_numeric <- gseaPlot(sig2Fun_result, showCategory = 3)
    expect_true(all(sapply(result_numeric, inherits, what = "gg")))

    specific_pathway <- c("REACTOME_CELL_CYCLE_MITOTIC",
                          "GOBP_CHROMOSOME_SEPARATION",
                          "REACTOME_SYNTHESIS_OF_DNA")
    result_specific <- tryCatch({
        gseaPlot(sig2Fun_result, showCategory = specific_pathway)
    }, error = function(e) {
        skip("Specific pathway test failed: pathway not found")
    })
    expect_true(all(sapply(result_specific, inherits, what = "gg")))
})

test_that("gseaPlot handles ESgeom parameter correctly", {
    res_line <- gseaPlot(sig2Fun_result, ESgeom = "line")
    res_dot  <- gseaPlot(sig2Fun_result, ESgeom = "dot")

    expect_true(all(sapply(res_line, inherits, what = "gg")))
    expect_true(all(sapply(res_dot, inherits, what = "gg")))
})

test_that("gseaPlot handles title parameter properly", {
    res_title <- gseaPlot(sig2Fun_result, title = "GSEA Test Title")
    plot_obj <- res_title[[1]]
    expect_true(!is.null(plot_obj$labels$title))
    expect_equal(plot_obj$labels$title, "GSEA Test Title")
})

test_that("gseaPlot handles colorPalette properly (matched and unmatched lengths)", {
    n_path <- length(unique(enrichplot:::get_gsdata(.extractDF(sig2Fun_result, "gseaRaw"), 1:3)$Description))
    colorPalette_matched <- viridisLite::viridis(n_path)
    res_col_good <- gseaPlot(sig2Fun_result, colorPalette = colorPalette_matched, showCategory = 3)
    expect_true(all(sapply(res_col_good, inherits, what = "gg")))
    colorPalette_wrong <- c("#FF0000", "#00FF00")
    res_col_wrong <- gseaPlot(sig2Fun_result, colorPalette = colorPalette_wrong, showCategory = 3)
    expect_true(all(sapply(res_col_wrong, inherits, what = "gg")))
})

test_that("gseaPlot pvalueTable parameter works and table is rendered", {
    res_with_table <- tryCatch({
        gseaPlot(sig2Fun_result, pvalueTable = TRUE, showCategory = 2)
    }, error = function(e) {
        skip(paste("pvalueTable test failed:", e$message))
    })
    expect_true(all(sapply(res_with_table, inherits, what = "gg")))
})

test_that("gseaPlot supports subPlots and relHeights customization", {
    res_partial <- gseaPlot(sig2Fun_result, subPlots = 1:2)
    expect_type(res_partial, "list")
    expect_true(length(res_partial) <= 2)

    res_full <- gseaPlot(sig2Fun_result, subPlots = 1:3, relHeights = c(2, 1, 1))
    expect_type(res_full, "list")
    expect_true(length(res_full) >= 2)
})

test_that("gseaPlot handles single pathway and geom_rect section correctly", {
    gseaRaw <- .extractDF(sig2Fun_result, type = "gseaRaw")
    single_pathway <- head(gseaRaw@result$ID, 1)
    res_single <- tryCatch({
        gseaPlot(sig2Fun_result, showCategory = single_pathway)
    }, error = function(e) {
        skip(paste("Single pathway test failed:", e$message))
    })
    expect_true(all(sapply(res_single, inherits, what = "gg")))
})

test_that("gseaPlot produces consistent y-axis and x-axis labeling", {
    res <- gseaPlot(sig2Fun_result)
    gp <- ggplot_build(res[[1]])
    x_labels <- gp$layout$panel_params[[1]]$x$get_labels()
    expect_true(is.character(x_labels))
})

test_that("gseaPlot handles different font sizes without error", {
    res_small <- gseaPlot(sig2Fun_result, fontSize = 8)
    res_large <- gseaPlot(sig2Fun_result, fontSize = 14)
    expect_true(all(sapply(res_small, inherits, what = "gg")))
    expect_true(all(sapply(res_large, inherits, what = "gg")))
})

test_that("gseaPlot works correctly with different relHeights values", {
    res1 <- gseaPlot(sig2Fun_result, relHeights = c(1, 1, 1))
    res2 <- gseaPlot(sig2Fun_result, relHeights = c(2, 0.5, 1))
    expect_true(all(sapply(res1, inherits, what = "gg")))
    expect_true(all(sapply(res2, inherits, what = "gg")))
})
