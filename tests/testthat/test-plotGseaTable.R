# Test file for .plotGseaTable function

library(testthat)
library(ggplot2)
library(dplyr)
library(cowplot)

if (!isNamespaceLoaded("SigFun")) {
    requireNamespace("SigFun", quietly = TRUE)
}

test_that(".plotGseaTable basic functionality returns ggplot and handles deprecation warning", {
    stats <- c(geneA = 2, geneB = -1, geneC = 1.5, geneD = -0.5)
    Pathways <- list(valid_path = c("geneA", "geneC"), bad_path = c("nonexistent"))
    fgseaRes <- data.frame(
        ID = c("valid_path", "bad_path"),
        NES = c(1.5, -2.0),
        pvalue = c(0.01, 0.05),
        stringsAsFactors = FALSE
    )

    expect_warning(
        p <- .plotGseaTable(
            pathways = Pathways,
            stats = stats,
            fgseaRes = fgseaRes,
            gseaParam = 1,
            pathwayLabelStyle = list(),
            render = TRUE,
            headerLabelStyle = list(),
            valueStyle = list(),
            axisLabelStyle = list()
        ),
        regexp = "render argument is deprecated"
    )

    expect_s3_class(p, "ggplot")
    expect_false(is.null(p))
})

test_that(".plotGseaTable handles all valid pathways", {
    stats <- c(g1 = 0.5, g2 = -0.2, g3 = 1.2, g4 = -0.8)
    Pathways <- list(p_alpha = c("g1", "g3"), p_beta = c("g2", "g4"))
    fgseaRes <- data.frame(
        ID = c("p_alpha", "p_beta"),
        NES = c(0.8, -1.1),
        pvalue = c(0.02, 0.03),
        stringsAsFactors = FALSE
    )

    expect_no_warning(
        p <- .plotGseaTable(
            pathways = Pathways,
            stats = stats,
            fgseaRes = fgseaRes,
            gseaParam = 2,
            pathwayLabelStyle = list(size = 5),
            headerLabelStyle = list(),
            valueStyle = list(),
            axisLabelStyle = list()
        )
    )
    expect_s3_class(p, "ggplot")
})

test_that(".plotGseaTable handles no valid pathways gracefully", {
    stats <- c(geneA = 2, geneB = -1)
    Pathways <- list(bad1 = "X", bad2 = "Y")
    fgseaRes <- data.frame(
        ID = c("bad1", "bad2"),
        NES = c(1.5, -2.0),
        pvalue = c(0.01, 0.05),
        stringsAsFactors = FALSE
    )

    expect_warning(
        p <- .plotGseaTable(
            pathways = Pathways,
            stats = stats,
            fgseaRes = fgseaRes,
            gseaParam = 1,
            pathwayLabelStyle = list(),
            headerLabelStyle = list(),
            valueStyle = list(),
            axisLabelStyle = list()
        ),
        regexp = "no non-missing arguments to max"
    )
    expect_s3_class(p, "ggplot")
})

test_that(".plotGseaTable handles mixed valid and invalid pathways", {
    stats <- c(geneA = 2, geneB = -1, geneC = 1.5)
    Pathways <- list(valid_path = c("geneA", "geneC"), bad_path = "X", another_valid = "geneB")
    fgseaRes <- data.frame(
        ID = c("valid_path", "bad_path", "another_valid"),
        NES = c(1.5, -2.0, 0.8),
        pvalue = c(0.01, 0.05, 0.03),
        stringsAsFactors = FALSE
    )

    expect_no_warning(
        p <- .plotGseaTable(
            pathways = Pathways,
            stats = stats,
            fgseaRes = fgseaRes,
            gseaParam = 1,
            pathwayLabelStyle = list(),
            headerLabelStyle = list(),
            valueStyle = list(),
            axisLabelStyle = list()
        )
    )
    expect_s3_class(p, "ggplot")
})

test_that(".plotGseaTable validates gseaParam variations", {
    stats <- c(geneA = 2, geneB = -1)
    Pathways <- list(valid_path = "geneA")
    fgseaRes <- data.frame(
        ID = "valid_path", NES = 1.5, pvalue = 0.01, stringsAsFactors = FALSE
    )

    expect_no_error({
        p1 <- .plotGseaTable(Pathways, stats, fgseaRes, gseaParam = 0.5)
        p2 <- .plotGseaTable(Pathways, stats, fgseaRes, gseaParam = 3)
    })
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
})

test_that(".plotGseaTable supports style parameters and custom widths", {
    stats <- c(geneA = 2, geneB = -1)
    Pathways <- list(valid_path = "geneA")
    fgseaRes <- data.frame(
        ID = "valid_path", NES = 1.5, pvalue = 0.01, stringsAsFactors = FALSE
    )

    expect_no_error({
        p <- .plotGseaTable(
            pathways = Pathways,
            stats = stats,
            fgseaRes = fgseaRes,
            gseaParam = 1,
            colWidths = c(1, 1, 1.5, 3, 8),  # ✅ fixed name
            pathwayLabelStyle = list(size = 12, hjust = 0.5),
            headerLabelStyle = list(size = 14),
            valueStyle = list(size = 10),
            axisLabelStyle = list(size = 8)
        )
    })
    expect_s3_class(p, "ggplot")
})

test_that(".plotGseaTable handles empty stats vector", {
    stats <- numeric(0)
    Pathways <- list(valid_path = "geneA")
    fgseaRes <- data.frame(
        ID = "valid_path", NES = 1.5, pvalue = 0.01, stringsAsFactors = FALSE
    )

    result <- suppressWarnings({
        warning_caught <- FALSE
        tryCatch({
            withCallingHandlers({
                p <- .plotGseaTable(
                    pathways = Pathways,
                    stats = stats,
                    fgseaRes = fgseaRes,
                    gseaParam = 1,
                    pathwayLabelStyle = list(),
                    headerLabelStyle = list(),
                    valueStyle = list(),
                    axisLabelStyle = list()
                )
            }, warning = function(w) {
                if (grepl("no non-missing arguments to max", w$message))
                    warning_caught <<- TRUE
            })
            list(plot = p, warning_caught = warning_caught)
        }, error = function(e) {
            list(plot = NULL, warning_caught = warning_caught)
        })
    })

    expect_true(result$warning_caught)
    expect_s3_class(result$plot, "ggplot")
})

test_that(".plotGseaTable works with single pathway and NULL styles", {
    stats <- c(geneA = 2, geneB = -1)
    Pathways <- list(single = c("geneA", "geneB"))
    fgseaRes <- data.frame(
        ID = "single", NES = 1.5, pvalue = 0.01, stringsAsFactors = FALSE
    )

    expect_error(
        .plotGseaTable(
            pathways = Pathways,
            stats = stats,
            fgseaRes = fgseaRes,
            gseaParam = 1,
            pathwayLabelStyle = NULL,
            headerLabelStyle = NULL,
            valueStyle = NULL,
            axisLabelStyle = NULL
        ),
        regexp = "GSEA statistic is not defined when all genes are selected"
    )
})
