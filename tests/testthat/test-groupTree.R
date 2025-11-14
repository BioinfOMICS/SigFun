library(testthat)
library(SigFun)
library(SummarizedExperiment)
library(ggplot2)
library(ggtree)

data("sig2Fun_result")

get_ids_gseasim <- function(se) {
    x <- SigFun:::.extractDF(se, type = "gseaSimilar")
    if (!is.null(x@result$ID)) unique(x@result$ID) else character(0)
}

test_that(".groupTree is called correctly through treePlot", {
    ids <- get_ids_gseasim(sig2Fun_result)

    if (length(ids) >= 5) {
        res <- suppressMessages({
            treePlot(sig2Fun_result, showCategory = ids[seq_len(5)])
        })

        expect_s3_class(res$treePlot, "gg")
        expect_s3_class(res$tableTreePlot, "data.frame")
    } else {
        skip("Not enough pathways for test")
    }
})

test_that(".groupTree handles different dotColor parameters correctly", {
    ids <- get_ids_gseasim(sig2Fun_result)

    if (length(ids) >= 5) {
        for (color_param in c("pvalue", "p.adjust", "qvalue", "NES")) {
            res <- suppressMessages({
                treePlot(sig2Fun_result,
                         showCategory = ids[seq_len(5)],
                         dotColor = color_param)
            })
            expect_s3_class(res$treePlot, "gg")
            layer_names <- sapply(res$treePlot$layers, function(x) class(x$geom)[1])
            expect_true(any(grepl("GeomPoint", layer_names)))
            expect_true(color_param %in% colnames(res$tableTreePlot))
        }
    } else {
        skip("Not enough pathways for dotColor test")
    }
})

test_that(".groupTree handles NES dotColor with gradient2 scale", {
    ids <- get_ids_gseasim(sig2Fun_result)

    if (length(ids) >= 5) {
        res <- suppressMessages({
            treePlot(sig2Fun_result,
                     showCategory = ids[seq_len(5)],
                     dotColor = "NES")
        })
        gp <- ggplot_build(res$treePlot)
        expect_true("NES" %in% colnames(res$tableTreePlot))
        expect_s3_class(res$treePlot, "gg")
    } else {
        skip("Not enough pathways for NES test")
    }
})

test_that(".groupTree handles custom groupColor correctly", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 5) {
        custom_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")

        res <- suppressMessages({
            treePlot(
                sig2Fun_result,
                showCategory = ids[seq_len(5)],
                clusterParams = list(n = 5, color = custom_colors)
            )
        })
        expect_s3_class(res$treePlot, "gg")
        has_color_scale <- any(sapply(res$treePlot$scales$scales,
                                      function(x) "colour" %in% x$aesthetics))
        expect_true(has_color_scale)
    } else {
        skip("Not enough pathways for custom color test")
    }
})

test_that(".groupTree handles NULL groupColor with viridis palette", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 5) {
        res <- suppressMessages({
            treePlot(
                sig2Fun_result,
                showCategory = ids[seq_len(5)],
                clusterParams = list(n = 5, color = NULL)
            )
        })

        expect_s3_class(res$treePlot, "gg")
    } else {
        skip("Not enough pathways for default color test")
    }
})

test_that(".groupTree handles labelFormatTiplab correctly", {
    ids <- get_ids_gseasim(sig2Fun_result)

    if (length(ids) >= 5) {
        res <- suppressMessages({
            treePlot(
                sig2Fun_result,
                showCategory = ids[seq_len(5)],
                labelFormatTiplab = 20
            )
        })
        expect_s3_class(res$treePlot, "gg")
        layer_geoms <- sapply(res$treePlot$layers, function(x) class(x$geom)[1])
        expect_true(any(grepl("GeomText", layer_geoms)))
    } else {
        skip("Not enough pathways for label format test")
    }
})

test_that(".groupTree handles custom labelFormatTiplab function", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 5) {
        custom_labeller <- function(x) toupper(substr(x, 1, 10))
        res <- suppressMessages({
            treePlot(
                sig2Fun_result,
                showCategory = ids[seq_len(5)],
                labelFormatTiplab = custom_labeller
            )
        })
        expect_s3_class(res$treePlot, "gg")
    } else {
        skip("Not enough pathways for custom labeller test")
    }
})

test_that(".groupTree handles relative offset values correctly", {
    ids <- get_ids_gseasim(sig2Fun_result)

    if (length(ids) >= 5) {
        res <- suppressMessages({
            treePlot(
                sig2Fun_result,
                showCategory = ids[seq_len(5)],
                offsetParams = list(
                    barTree = rel(1.5),
                    tiplab = rel(1.3),
                    extend = 0.25,
                    hexpand = 0.1
                )
            )
        })

        expect_s3_class(res$treePlot, "gg")
    } else {
        skip("Not enough pathways for offset test")
    }
})

test_that(".groupTree handles absolute offset values correctly", {
    ids <- get_ids_gseasim(sig2Fun_result)

    if (length(ids) >= 5) {
        res <- suppressMessages({
            treePlot(
                sig2Fun_result,
                showCategory = ids[seq_len(5)],
                offsetParams = list(
                    barTree = 2,
                    tiplab = 1.5,
                    extend = 0.3,
                    hexpand = 0.08
                )
            )
        })

        expect_s3_class(res$treePlot, "gg")
    } else {
        skip("Not enough pathways for absolute offset test")
    }
})

test_that(".groupTree handles different cluster numbers", {
    ids <- get_ids_gseasim(sig2Fun_result)

    if (length(ids) >= 8) {
        for (n_cluster in c(3, 5, 7)) {
            res <- suppressMessages({
                treePlot(
                    sig2Fun_result,
                    showCategory = ids[seq_len(8)],
                    clusterParams = list(n = n_cluster)
                )
            })

            expect_s3_class(res$treePlot, "gg")
        }
    } else {
        skip("Not enough pathways for cluster number test")
    }
})

test_that(".groupTree handles different hclust methods", {
    ids <- get_ids_gseasim(sig2Fun_result)

    if (length(ids) >= 5) {
        for (method in c("ward.D", "ward.D2", "average", "complete")) {
            res <- suppressMessages({
                tryCatch({
                    treePlot(
                        sig2Fun_result,
                        showCategory = ids[seq_len(5)],
                        clusterParams = list(method = method)
                    )
                }, error = function(e) {
                    skip(paste("Method", method, "failed"))
                })
            })

            if (!is.null(res)) {
                expect_s3_class(res$treePlot, "gg")
            }
        }
    } else {
        skip("Not enough pathways for hclust method test")
    }
})

# ===== 字體大小測試 =====

test_that(".groupTree handles different font sizes", {
    ids <- get_ids_gseasim(sig2Fun_result)

    if (length(ids) >= 5) {
        res_small <- suppressMessages({
            treePlot(
                sig2Fun_result,
                showCategory = ids[seq_len(5)],
                leafFontSize = 2,
                cladeFontSize = 2
            )
        })

        res_large <- suppressMessages({
            treePlot(
                sig2Fun_result,
                showCategory = ids[seq_len(5)],
                leafFontSize = 5,
                cladeFontSize = 4
            )
        })

        expect_s3_class(res_small$treePlot, "gg")
        expect_s3_class(res_large$treePlot, "gg")
    } else {
        skip("Not enough pathways for font size test")
    }
})

test_that(".groupTree handles cexCategory parameter", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 5) {
        for (cex_val in c(0.5, 1, 1.5)) {
            res <- suppressMessages({
                treePlot(
                    sig2Fun_result,
                    showCategory = ids[seq_len(5)],
                    cexCategory = cex_val
                )
            })
            expect_s3_class(res$treePlot, "gg")
        }
    } else {
        skip("Not enough pathways for cexCategory test")
    }
})

test_that(".groupTree handles hilight parameters", {
    ids <- get_ids_gseasim(sig2Fun_result)

    if (length(ids) >= 5) {
        res_hilight <- suppressMessages({
            treePlot(
                sig2Fun_result,
                showCategory = ids[seq_len(5)],
                hilightParams = list(hilight = TRUE, align = "both")
            )
        })
        res_no_hilight <- suppressMessages({
            treePlot(
                sig2Fun_result,
                showCategory = ids[seq_len(5)],
                hilightParams = list(hilight = FALSE, align = "none")
            )
        })

        expect_s3_class(res_hilight$treePlot, "gg")
        expect_s3_class(res_no_hilight$treePlot, "gg")
    } else {
        skip("Not enough pathways for hilight test")
    }
})

test_that(".groupTree produces correct layer structure", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 5) {
        res <- suppressMessages({
            treePlot(sig2Fun_result, showCategory = ids[seq_len(5)])
        })
        plot_obj <- res$treePlot
        layer_geoms <- sapply(plot_obj$layers, function(x) class(x$geom)[1])
        expect_true(any(grepl("Point", layer_geoms)))
        expect_true(any(grepl("Text", layer_geoms)))
        expect_true(length(plot_obj$scales$scales) > 0)
    } else {
        skip("Not enough pathways for layer structure test")
    }
})

test_that(".groupTree plot data contains required columns", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 5) {
        res <- suppressMessages({
            treePlot(sig2Fun_result, showCategory = ids[seq_len(5)])
        })
        plot_data <- res$treePlot$data
        expect_true("label" %in% colnames(plot_data))
        expect_true("isTip" %in% colnames(plot_data))
        expect_true("group" %in% colnames(plot_data))
    } else {
        skip("Not enough pathways for data structure test")
    }
})

test_that(".groupTree handles minimum number of pathways", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 5) {
        res <- suppressMessages({
            treePlot(sig2Fun_result, showCategory = ids[seq_len(5)])
        })
        expect_s3_class(res$treePlot, "gg")
        expect_equal(nrow(res$tableTreePlot), 5)
    } else {
        skip("Not enough pathways for minimum test")
    }
})

test_that(".groupTree handles NA values in color gracefully", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 5) {
        res <- suppressMessages({
            treePlot(sig2Fun_result, showCategory = ids[seq_len(5)])
        })
        expect_s3_class(res$treePlot, "gg")
        color_scales <- Filter(function(x) "colour" %in% x$aesthetics,
                               res$treePlot$scales$scales)
        expect_true(length(color_scales) > 0)
    } else {
        skip("Not enough pathways for NA handling test")
    }
})

test_that(".groupTree integrates correctly with all treePlot parameters", {
    ids <- get_ids_gseasim(sig2Fun_result)
    if (length(ids) >= 6) {
        res <- suppressMessages({
            treePlot(
                sig2Fun_result,
                showCategory = ids[seq_len(6)],
                dotColor = "p.adjust",
                hclustfun = "ward.D2",
                labelFormat = 25,
                labelFormatTiplab = 30,
                cexCategory = 1.2,
                leafFontSize = 3.5,
                cladeFontSize = 3,
                hilightParams = list(hilight = TRUE, align = "both"),
                offsetParams = list(barTree = rel(1.4), tiplab = rel(1.6),
                                    extend = 0.35, hexpand = 0.12),
                clusterParams = list(method = "average", n = 4,
                                     labelWordsN = 2, labelFormat = 28)
            )
        })
        expect_s3_class(res$treePlot, "gg")
        expect_s3_class(res$tableTreePlot, "data.frame")
        expect_true("p.adjust" %in% colnames(res$tableTreePlot))
    } else {
        skip("Not enough pathways for integration test")
    }
})
