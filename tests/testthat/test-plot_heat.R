#' Unit Tests for GSEA Plot Functions
#'
#' @author Claude
#' @date April 24, 2025

library(testthat)
library(SummarizedExperiment)
library(magrittr)
library(dplyr)
library(S4Vectors)
library(ggplot2)
library(data.table)
library(tibble)
library(cowplot)
library(fgsea)

# ------------------------------------------------------------------------------
# Test fixtures and helper functions
# ------------------------------------------------------------------------------

# Create mock data for tests
create_mock_data <- function() {
    # Create mock fgsea result
    mock_fgsea_result <- data.frame(
        pathway = c(
            "KEGG_pathway1", "KEGG_pathway2", "GOBP_pathway1",
            "REACTOME_pathway1", "WP_pathway1", "HALLMARK_pathway1"
        ),
        pval = c(0.01, 0.02, 0.03, 0.04, 0.001, 0.02),
        qval = c(0.02, 0.03, 0.04, 0.05, 0.01, 0.03),
        NES = c(2.5, 2.0, 1.5, 1.2, 3.0, -2.2),
        leadingEdge = c(
            "Gene1|Gene2", "Gene3|Gene4", "Gene5|Gene6",
            "Gene7|Gene8", "Gene9|Gene10", "Gene11|Gene12"
        )
    )

    # Create mock correlation data
    mock_cor_df <- data.frame(
        gene = paste0("ENSG", sprintf("%011d", 1:20)),
        cor = seq(-0.9, 0.9, length.out = 20)
    )

    # Create mock row data for SummarizedExperiment
    mock_row_data <- data.frame(
        gene_symbol = paste0("Gene", 1:20),
        ensg_id = paste0("ENSG", sprintf("%011d", 1:20)),
        gene_biotype = rep(c("protein_coding", "lncRNA"), each = 10)
    )

    # Create mock pathways
    mock_pathways <- list(
        "KEGG_pathway1" = c("Gene1", "Gene2", "Gene3"),
        "KEGG_pathway2" = c("Gene4", "Gene5", "Gene6"),
        "GOBP_pathway1" = c("Gene7", "Gene8", "Gene9"),
        "REACTOME_pathway1" = c("Gene10", "Gene11", "Gene12"),
        "WP_pathway1" = c("Gene13", "Gene14", "Gene15"),
        "HALLMARK_pathway1" = c("Gene16", "Gene17", "Gene18")
    )

    # Create mock stats
    mock_stats <- setNames(
        seq(-0.9, 0.9, length.out = 20),
        paste0("Gene", 1:20)
    )

    # Create a mock SummarizedExperiment object
    mock_SE <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:100, nrow = 20, ncol = 5)),
        rowData = mock_row_data
    )
    S4Vectors::metadata(mock_SE)$fgsea <- mock_fgsea_result
    S4Vectors::metadata(mock_SE)$cor.df <- mock_cor_df

    return(list(
        SE = mock_SE,
        pathways = mock_pathways,
        stats = mock_stats,
        fgsea_result = mock_fgsea_result
    ))
}

# ------------------------------------------------------------------------------
# Tests for .cowplotText function
# ------------------------------------------------------------------------------

test_that(".cowplotText creates a ggplot object with text", {
    # Test with default style
    p1 <- .cowplotText("Test Text", list(size = 12))
    expect_s3_class(p1, "ggplot")

    # Test with custom style
    custom_style <- list(size = 14, color = "red", fontface = "bold")
    p2 <- .cowplotText("Custom Text", custom_style)
    expect_s3_class(p2, "ggplot")
})

# ------------------------------------------------------------------------------
# Tests for .plotEnrichmentData function
# ------------------------------------------------------------------------------

test_that(".plotEnrichmentData calculates enrichment data correctly", {
    # Create test data
    mock_data <- create_mock_data()
    pathway <- c("Gene1", "Gene5", "Gene10")
    stats <- mock_data$stats

    # Run function
    result <- .plotEnrichmentData(pathway, stats)

    # Check structure of result
    expect_type(result, "list")
    expect_named(result, c("curve", "ticks", "stats", "posES", "negES", "spreadES", "maxAbsStat"))

    # Check data types
    expect_s3_class(result$curve, "data.table")
    expect_s3_class(result$ticks, "data.table")
    expect_s3_class(result$stats, "data.table")
    expect_type(result$posES, "double")
    expect_type(result$negES, "double")
    expect_type(result$spreadES, "double")
    expect_type(result$maxAbsStat, "double")

    # Check for errors with invalid inputs
    expect_error(.plotEnrichmentData(pathway, c(stats, Inf)), "Not all stats values are finite numbers")
})

# ------------------------------------------------------------------------------
# Tests for .plotEnrichment function
# ------------------------------------------------------------------------------

test_that(".plotEnrichment creates a ggplot visualization", {
    # Create test data
    mock_data <- create_mock_data()
    pathway <- c("Gene1", "Gene5", "Gene10")
    stats <- mock_data$stats

    # Run function
    plot <- .plotEnrichment(pathway, stats, color = 0.5)

    # Check result type
    expect_s3_class(plot, "ggplot")

    # Test with negative color
    plot_neg <- .plotEnrichment(pathway, stats, color = -0.5)
    expect_s3_class(plot_neg, "ggplot")

    # Check for errors with invalid inputs
    expect_error(.plotEnrichment(pathway, c(stats, Inf), color = 0.5))
})

# ------------------------------------------------------------------------------
# Tests for .get_ps function
# ------------------------------------------------------------------------------

test_that(".get_ps creates a list of plots", {
    skip_if_not_installed("cowplot")

    # Create test data
    mock_data <- create_mock_data()

    # Select a subset of pathways for testing
    test_pathways <- mock_data$pathways[1:2]

    # Run function
    results <- .get_ps(
        Pathways = test_pathways,
        fgseaRes = mock_data$fgsea_result,
        maxNES = 3.0,
        valueStyle = list(size = 12),
        headerLabelStyle = list(size = 12),
        pathways.all = mock_data$pathways,
        pathwayLabelStyleDefault = list(size = 12, hjust = 0, x = 0.05, vjust = 0.5),
        stats = mock_data$stats
    )

    # Check structure
    expect_type(results, "list")
    expect_equal(length(results), length(test_pathways))

    # Check content of each pathway result
    for (i in seq_along(results)) {
        pathway_results <- results[[i]]
        expect_type(pathway_results, "list")
        expect_equal(length(pathway_results), 5)  # NES bar, NES value, p-value, enrichment plot, pathway name

        # Check each plot in the results
        expect_s3_class(pathway_results[[1]], "ggplot")  # NES bar
        expect_s3_class(pathway_results[[2]], "ggplot")  # NES value text
        expect_s3_class(pathway_results[[3]], "ggplot")  # p-value text
        expect_s3_class(pathway_results[[4]], "ggplot")  # enrichment plot
        expect_s3_class(pathway_results[[5]], "ggplot")  # pathway name text
    }
})

# ------------------------------------------------------------------------------
# Tests for .plotGseaTable function
# ------------------------------------------------------------------------------

test_that(".plotGseaTable creates a plot grid", {
    skip_if_not_installed("cowplot")

    # Create test data
    mock_data <- create_mock_data()

    # Run function
    plot <- .plotGseaTable(
        pathways.all = mock_data$pathways,
        Pathways = mock_data$pathways[1:2],
        stats = mock_data$stats,
        fgseaRes = mock_data$fgsea_result
    )

    # Check result type
    expect_s3_class(plot, "ggplot")

    # Test with custom styling
    plot_custom <- .plotGseaTable(
        pathways.all = mock_data$pathways,
        Pathways = mock_data$pathways[1:2],
        stats = mock_data$stats,
        fgseaRes = mock_data$fgsea_result,
        pathwayLabelStyle = list(size = 14, color = "blue"),
        headerLabelStyle = list(size = 14, color = "red"),
        valueStyle = list(size = 14, color = "green")
    )
    expect_s3_class(plot_custom, "ggplot")

    # Test with deprecated render parameter
    expect_warning(
        .plotGseaTable(
            pathways.all = mock_data$pathways,
            Pathways = mock_data$pathways[1:2],
            stats = mock_data$stats,
            fgseaRes = mock_data$fgsea_result,
            render = TRUE
        ),
        "render argument is deprecated"
    )

    # Test with empty pathway (should not error)
    suppressWarnings(plot_empty <- .plotGseaTable(
        pathways.all = mock_data$pathways,
        Pathways = list("empty_pathway" = character(0)),
        stats = mock_data$stats,
        fgseaRes = mock_data$fgsea_result
    ))
    expect_s3_class(plot_empty, "ggplot")
})

# ------------------------------------------------------------------------------
# Tests for plot_heat function
# ------------------------------------------------------------------------------

test_that("plot_heat processes data and returns expected output", {
    skip_if_not_installed("cowplot")

    # Create mock data
    mock_data <- create_mock_data()

    # Set up a mock for .plotGseaTable function since it's complex
    mock_plotGseaTable <- function(pathways.all, Pathways, stats, fgseaRes, ...) {
        ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::annotate(geom = "text", x = 1, y = 1, label = "Mock GSEA table")
    }

    # Temporarily assign the mock function
    orig_plotGseaTable <- .plotGseaTable
    assign(".plotGseaTable", mock_plotGseaTable, envir = .GlobalEnv)

    # Run function
    temp_dir <- tempdir()
    data("demo_GSE181574")
    SE_data.cor <- sigCor(SE_data=SE_GSE181574, cor.method="logit",
                          output_path=temp_dir, Z.transform=FALSE)
    metadata(SE_data.cor) <- list(cor.df=metadata(SE_data.cor)$cor.df)
    suppressWarnings(SE_data.fgsea <- sig2GSEA(SE_data.cor=SE_data.cor,
                              ranking.method="stat", output_path=temp_dir,
                              pathways.all=pathways.all))
    metadata(SE_data.fgsea) <- list(fgsea=metadata(SE_data.fgsea)$fgsea,
                                    cor.df=metadata(SE_data.cor)$cor.df)
    result <- plot_heat(
        SE_data.fgsea = SE_data.fgsea,
        output_path = temp_dir,
        significat_type = "pval",
        strings = c("KEGG", "GOBP", "REACTOME", "WP", "HALLMARK"),
        topN = 5,
        pathways.all = pathways.all,
        ranking.method = "stat",
        plot_out = FALSE
    )

    # Restore the original function
    assign(".plotGseaTable", orig_plotGseaTable, envir = .GlobalEnv)

    # Check result structure
    expect_type(result, "list")
    expect_equal(length(result), 5)  # Should match number of strings
    expect_true(all(names(result) %in% c("KEGG", "GOBP", "REACTOME", "WP", "HALLMARK")))

    # Check each returned plot is a ggplot object
    for (plot_name in names(result)) {
        expect_s3_class(result[[plot_name]], "ggplot")
    }
})

test_that("plot_heat handles different significance types", {
    skip_if_not_installed("cowplot")

    # Create mock data
    mock_data <- create_mock_data()

    # Mock .plotGseaTable
    mock_plotGseaTable <- function(pathways.all, Pathways, stats, fgseaRes, ...) {
        ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::annotate(geom = "text", x = 1, y = 1, label = "Mock GSEA table")
    }

    # Temporarily assign the mock function
    orig_plotGseaTable <- .plotGseaTable
    assign(".plotGseaTable", mock_plotGseaTable, envir = .GlobalEnv)

    # Run with p-value
    temp_dir <- tempdir()
    p_result <- plot_heat(
        SE_data.fgsea = mock_data$SE,
        output_path = temp_dir,
        significat_type = "pval",
        strings = c("KEGG"),
        topN = 10,
        pathways.all = mock_data$pathways,
        ranking.method = "stat",
        plot_out = FALSE
    )

    # Run with q-value
    q_result <- plot_heat(
        SE_data.fgsea = mock_data$SE,
        output_path = temp_dir,
        significat_type = "qval",
        strings = c("KEGG"),
        topN = 10,
        pathways.all = mock_data$pathways,
        ranking.method = "stat",
        plot_out = FALSE
    )

    # Restore the original function
    assign(".plotGseaTable", orig_plotGseaTable, envir = .GlobalEnv)

    # Both should return ggplot objects
    expect_s3_class(p_result[["KEGG"]], "ggplot")
    expect_s3_class(q_result[["KEGG"]], "ggplot")
})

test_that("plot_heat handles different ranking methods", {
    skip_if_not_installed("cowplot")

    # Create mock data
    mock_data <- create_mock_data()

    # Mock .plotGseaTable
    mock_plotGseaTable <- function(pathways.all, Pathways, stats, fgseaRes, ...) {
        ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::annotate(geom = "text", x = 1, y = 1, label = "Mock GSEA table")
    }

    # Temporarily assign the mock function
    orig_plotGseaTable <- .plotGseaTable
    assign(".plotGseaTable", mock_plotGseaTable, envir = .GlobalEnv)

    # Run with stat ranking
    temp_dir <- tempdir()
    stat_result <- plot_heat(
        SE_data.fgsea = mock_data$SE,
        output_path = temp_dir,
        strings = c("KEGG"),
        topN = 10,
        pathways.all = mock_data$pathways,
        ranking.method = "stat",
        plot_out = FALSE
    )

    # Run with absolute stat ranking
    abs_stat_result <- plot_heat(
        SE_data.fgsea = mock_data$SE,
        output_path = temp_dir,
        strings = c("KEGG"),
        topN = 10,
        pathways.all = mock_data$pathways,
        ranking.method = "abs.stat",
        plot_out = FALSE
    )

    # Restore the original function
    assign(".plotGseaTable", orig_plotGseaTable, envir = .GlobalEnv)

    # Both should return ggplot objects
    expect_s3_class(stat_result[["KEGG"]], "ggplot")
    expect_s3_class(abs_stat_result[["KEGG"]], "ggplot")
})

test_that("plot_heat throws appropriate errors for invalid inputs", {
    # Create a mock SummarizedExperiment object WITHOUT fgsea metadata
    mock_SE_no_fgsea <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:100, nrow = 20, ncol = 5)),
        rowData = data.frame(
            gene_symbol = paste0("Gene", 1:20),
            ensg_id = paste0("ENSG", sprintf("%011d", 1:20)),
            gene_biotype = rep(c("protein_coding", "lncRNA"), each = 10)
        )
    )

    # Set up temporary directory
    temp_dir <- tempdir()

    # Test with missing fgsea metadata
    expect_error(
        plot_heat(
            SE_data.fgsea = mock_SE_no_fgsea,
            output_path = temp_dir,
            pathways.all = list()
        )
    )

    # Create a mock with fgsea but no cor.df
    mock_SE_no_cor <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:100, nrow = 20, ncol = 5)),
        rowData = data.frame(
            gene_symbol = paste0("Gene", 1:20),
            ensg_id = paste0("ENSG", sprintf("%011d", 1:20)),
            gene_biotype = rep(c("protein_coding", "lncRNA"), each = 10)
        )
    )

    mock_fgsea_result <- data.frame(
        pathway = c("KEGG_pathway1"),
        pval = 0.01,
        qval = 0.02,
        NES = 2.5,
        leadingEdge = "Gene1|Gene2"
    )

    S4Vectors::metadata(mock_SE_no_cor)$fgsea <- mock_fgsea_result

    # Test with missing cor.df metadata
    expect_error(
        plot_heat(
            SE_data.fgsea = mock_SE_no_cor,
            output_path = temp_dir,
            pathways.all = list()
        )
    )

    # Create a mock SummarizedExperiment without necessary columns in rowData
    mock_row_data_incomplete <- data.frame(
        name = paste0("Gene", 1:20)  # Missing required columns
    )

    mock_SE_bad_rowdata <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:100, nrow = 20, ncol = 5)),
        rowData = mock_row_data_incomplete
    )

    S4Vectors::metadata(mock_SE_bad_rowdata)$fgsea <- mock_fgsea_result
    S4Vectors::metadata(mock_SE_bad_rowdata)$cor.df <- data.frame(
        gene = paste0("ENSG", sprintf("%011d", 1:20)),
        cor = seq(-0.9, 0.9, length.out = 20)
    )

    # Test with incomplete rowData
    expect_error(
        plot_heat(
            SE_data.fgsea = mock_SE_bad_rowdata,
            output_path = temp_dir,
            pathways.all = list()
        )
    )
})
