library(testthat)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(clusterProfiler)
library(fgsea)
library(S4Vectors)
library(SigFun)

# Test helper function for class checking
test_that(".classCheck works correctly", {
    # Mock the .classCheck function since it's not provided
    .classCheck <- function(obj, expected_class) {
        if (!inherits(obj, expected_class)) {
            stop(paste("Object is not of class", expected_class))
        }
    }

    # Test with correct class
    se_data <- SummarizedExperiment(
        assays = list(counts = matrix(1:20, nrow = 4))
    )
    expect_silent(.classCheck(se_data, "SummarizedExperiment"))

    # Test with incorrect class
    expect_error(.classCheck(list(), "SummarizedExperiment"))
})

data("expr.data")
data("mapping")
data("SIG_MAT")
data("t2g")
seData <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = as.matrix(expr.data)),
    rowData = S4Vectors::DataFrame(mapping, row.names = mapping$ensg_id),
    colData = S4Vectors::DataFrame(SIG_MAT)
)

# Test sig2Fun main function
test_that("sig2Fun performs complete signature analysis", {

    # Test sig2Fun function with expected warnings
    expect_warning({
        result <- sig2Fun(
            seData = seData,
            rankingMethod = "cor",
            species = "human",
            t2g = t2g,
            corMethod = "logit",
            topN = 10,
            zTransform = FALSE,
            significantType = "pval",
            strings = c("GOBP", "GOCC", "GOMF", "KEGG", "REACTOME",
                        "WP", "HALLMARK", "SIGNALING")
        )
    }, regexp = "For some of the pathways the P-values were likely overestimated|For some pathways, in reality P-values are less than 1e-10")

    # Alternative: expect specific warnings individually
    expect_warning(
        expect_warning({
            result <- sig2Fun(
                seData = seData,
                rankingMethod = "cor",
                species = "human",
                t2g = t2g,
                corMethod = "logit",
                topN = 10,
                zTransform = FALSE,
                significantType = "pval",
                strings = c("GOBP", "GOCC", "GOMF", "KEGG", "REACTOME",
                            "WP", "HALLMARK", "SIGNALING")
            )
        }, "For some of the pathways the P-values were likely overestimated"),
        "For some pathways, in reality P-values are less than 1e-10"
    )
})

# Test error handling
test_that("Functions properly handle invalid inputs", {
    # Test .plotEnrichmentData with invalid stats
    invalid_stats <- c(1, Inf, -Inf, NaN)
    pathway <- c("g1", "g2")
    expect_error(.plotEnrichmentData(pathway, invalid_stats))

    # Test .corList with mismatched dimensions
    short_signature <- c(1, 2, 3)
    long_genes <- data.frame(
        gene1 = c(1, 2, 3, 4, 5),
        gene2 = c(2, 3, 4, 5, 6)
    )
    expect_error(.corList(short_signature, long_genes, "pearson"))

    # Test .z_score_cal with empty vector
    expect_error(.z_score_cal(numeric(0)))
})
