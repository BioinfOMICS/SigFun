library(testthat)
library(SigFun)
library(dplyr)
library(SummarizedExperiment)

data("sig2Fun_result")
data("t2g")
pathwaysAll <- split(t2g$ensembl_gene, t2g$gs_name)

test_that("plot_heat basic functionality works correctly", {

    result <- plot_heat(
        seDataFgsea = sig2Fun_result,
        pathwaysAll = pathwaysAll
    )
    expect_type(result, "list")
    expected_names <- c("GOBP", "GOCC", "GOMF", "KEGG", "REACTOME", "WP", "HALLMARK")
    expect_true(all(expected_names %in% names(result)))
    expect_true(all(sapply(result, function(x) inherits(x, "ggplot") || is.list(x))))
})


test_that("plot_heat handles custom parameters properly", {
    result_custom <- plot_heat(
        seDataFgsea = sig2Fun_result,
        pathwaysAll = pathwaysAll,
        significantType = "padj",
        topN = 5,
        rankingMethod = "cor",
        colWidths = c(1.5, 1, 1, 5, 4),
        fontSize = 9,
        pvalueType = "p.adjust"
    )
    expect_type(result_custom, "list")
    expect_named(result_custom)
})


test_that("plot_heat throws informative errors for invalid parameters", {
    expect_error(
        plot_heat(sig2Fun_result, pathwaysAll, colWidths = c(1, 2)),
        regexp = "must be a vector of length 5"
    )
    expect_error(
        plot_heat(sig2Fun_result, pathwaysAll, rankingMethod = "NotExistCol"),
        regexp = "must be included in"
    )
    expect_error(
        plot_heat(sig2Fun_result, pathwaysAll, pvalueType = "invalid"),
        regexp = "'arg' should be one of"
    )
})
