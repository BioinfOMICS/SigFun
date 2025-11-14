library(testthat)
library(SigFun)
library(SummarizedExperiment)

data("sig2Fun_result")

get_ids_gseasim <- function(se) {
    x <- SigFun:::.extractDF(se, type = "gseaSimilar")
    if (!is.null(x@result$ID)) unique(x@result$ID) else character(0)
}

test_that("treePlot basic output structure", {
    res <- treePlot(sig2Fun_result)
    expect_type(res, "list")
    expect_named(res, c("treePlot", "tableTreePlot"))
    expect_s3_class(res$treePlot, "gg")
    expect_s3_class(res$tableTreePlot, "data.frame")
    expect_true(all(c("labelName", "count") %in% colnames(res$tableTreePlot)))
})

test_that("treePlot works with numeric showCategory >= 5", {
    res <- treePlot(sig2Fun_result, showCategory = 6)
    expect_s3_class(res$treePlot, "gg")
    expect_lte(nrow(res$tableTreePlot), 6)
})

test_that("treePlot works with specific pathway vector (>= 5)", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 5) {
        pick <- ids[seq_len(5)]
        res <- treePlot(sig2Fun_result, showCategory = pick)
        expect_s3_class(res$treePlot, "gg")
        expect_lte(nrow(res$tableTreePlot), length(pick))
    } else {
        skip("Not enough pathways in gseaSimilar@result$ID")
    }
})

test_that("treePlot supports dotColor options", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 5) {
        for (colv in c("pvalue", "p.adjust", "qvalue", "NES")) {
            res <- treePlot(sig2Fun_result, showCategory = ids[seq_len(5)], dotColor = colv)
            expect_s3_class(res$treePlot, "gg")
            expect_true(colv %in% colnames(res$tableTreePlot))
        }
    } else {
        skip("Not enough pathways for dotColor test")
    }
})

test_that("treePlot accepts clusterParams override (method, n, color, labelWordsN, labelFormat)", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 5) {
        pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
        res <- treePlot(
            sig2Fun_result,
            showCategory = ids[seq_len(5)],
            clusterParams = list(method = "average", n = 5, color = pal,
                                 labelWordsN = 3, labelFormat = 20)
        )
        expect_s3_class(res$treePlot, "gg")
    } else {
        skip("Not enough pathways for clusterParams test")
    }
})

test_that("treePlot accepts hilightParams / offsetParams", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 5) {
        res <- treePlot(
            sig2Fun_result,
            showCategory = ids[seq_len(5)],
            hilightParams = list(hilight = TRUE, align = "both"),
            offsetParams  = list(barTree = rel(1.2), tiplab = rel(1.4),
                                 extend = 0.25, hexpand = 0.08)
        )
        expect_s3_class(res$treePlot, "gg")
    } else {
        skip("Not enough pathways for hilight/offset test")
    }
})

test_that("treePlot respects labelFormat / labelFormatTiplab", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 5) {
        res <- treePlot(
            sig2Fun_result,
            showCategory = ids[seq_len(5)],
            labelFormat = 25,
            labelFormatTiplab = 30
        )
        expect_s3_class(res$treePlot, "gg")
    } else {
        skip("Not enough pathways for labelFormat test")
    }
})

test_that("treePlot error when numeric showCategory < 5", {
    expect_error(
        treePlot(sig2Fun_result, showCategory = 2),
        regexp = "must be larger than"
    )
})

test_that("treePlot error when pathway vector length < 5 or not found", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 2) {
        expect_error(
            treePlot(sig2Fun_result, showCategory = ids[seq_len(2)]),
            regexp = "must provide at least"
        )
    } else {
        skip("Not enough pathways to test <5 vector error")
    }
})

test_that("treePlot table columns by dotColor are present", {
    res <- treePlot(sig2Fun_result)
    expect_true(any(c("pvalue", "p.adjust", "qvalue", "NES") %in% colnames(res$tableTreePlot)))
})

test_that("treePlot produces valid sizes and non-negative counts", {
    res <- treePlot(sig2Fun_result)
    expect_true(all(is.finite(res$tableTreePlot$count)))
    expect_true(all(res$tableTreePlot$count >= 0))
})
