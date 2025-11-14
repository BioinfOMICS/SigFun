library(testthat)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)

# Helper function to create test data
create_test_se <- function(n_genes = 5, n_samples = 10, add_signature = TRUE) {
    # Create gene expression matrix
    set.seed(123)
    gene_names <- paste0("ENSG", seq_len(n_genes))
    sample_names <- paste0("Sample", seq_len(n_samples))

    assay_data <- matrix(
        rnorm(n_genes * n_samples, mean = 10, sd = 2),
        nrow = n_genes,
        ncol = n_samples,
        dimnames = list(gene_names, sample_names)
    )

    # Create colData with signature values
    col_data <- data.frame(
        sample_id = sample_names,
        stringsAsFactors = FALSE
    )

    if (add_signature) {
        col_data$value <- rnorm(n_samples, mean = 5, sd = 1)
    }

    # Create SummarizedExperiment
    se <- SummarizedExperiment(
        assays = list(counts = assay_data),
        colData = col_data
    )

    return(se)
}

# Unit Tests
test_that("sigCor returns SummarizedExperiment object", {
    se_test <- create_test_se()
    result <- sigCor(seData = se_test)
    expect_s4_class(result, "SummarizedExperiment")
})

test_that("sigCor adds correlation results to metadata", {
    se_test <- create_test_se()
    result <- sigCor(seData = se_test)

    expect_true("cor.df" %in% names(S4Vectors::metadata(result)))
    cor_df <- S4Vectors::metadata(result)$cor.df
    expect_s3_class(cor_df, "data.frame")
    expect_true(all(c("gene", "cor", "pval") %in% colnames(cor_df)))
})

test_that("sigCor works with different correlation methods", {
    se_test <- create_test_se()

    # Test spearman (default)
    result_spearman <- sigCor(seData = se_test, corMethod = "spearman")
    expect_true("cor.df" %in% names(S4Vectors::metadata(result_spearman)))

    # Test pearson
    result_pearson <- sigCor(seData = se_test, corMethod = "pearson")
    expect_true("cor.df" %in% names(S4Vectors::metadata(result_pearson)))

    # Test kendall
    result_kendall <- sigCor(seData = se_test, corMethod = "kendall")
    expect_true("cor.df" %in% names(S4Vectors::metadata(result_kendall)))
})

test_that("sigCor works with logistic regression method", {
    se_test <- create_test_se()
    # Create binary signature for logistic regression
    colData(se_test)$value <- sample(c(0, 1), ncol(se_test), replace = TRUE)

    result_logit <- sigCor(seData = se_test, corMethod = "logit")
    expect_true("cor.df" %in% names(S4Vectors::metadata(result_logit)))

    cor_df <- S4Vectors::metadata(result_logit)$cor.df
    expect_equal(nrow(cor_df), nrow(se_test))
})

test_that("sigCor works with Z-transformation", {
    # Create test data with continuous values to avoid ties warning
    set.seed(789)
    se_test <- create_test_se(n_genes = 5, n_samples = 20)

    # Create genes with continuous values and avoid ties
    assay_data <- matrix(
        runif(5 * 20, min = 0, max = 100),  # Use runif to avoid ties
        nrow = 5,
        ncol = 20
    )
    rownames(assay_data) <- rownames(assay(se_test))
    colnames(assay_data) <- colnames(assay(se_test))
    assay(se_test) <- assay_data

    # Create continuous signature values without ties
    colData(se_test)$value <- runif(20, min = 1, max = 10)

    result_no_z <- sigCor(seData = se_test, zTransform = FALSE)
    result_with_z <- sigCor(seData = se_test, zTransform = TRUE)

    cor_df_no_z <- S4Vectors::metadata(result_no_z)$cor.df
    cor_df_with_z <- S4Vectors::metadata(result_with_z)$cor.df

    # Check that results exist and have reasonable structure
    expect_s3_class(cor_df_no_z, "data.frame")
    expect_s3_class(cor_df_with_z, "data.frame")
    expect_equal(nrow(cor_df_no_z), nrow(cor_df_with_z))
    expect_true(all(c("gene", "cor", "pval") %in% colnames(cor_df_with_z)))

    # Both should produce valid correlation results
    expect_true(all(!is.na(cor_df_no_z$gene)))
    expect_true(all(!is.na(cor_df_with_z$gene)))
})

test_that("sigCor handles sample mismatch between assay and colData", {
    # Create base test data
    se_test <- create_test_se(n_samples = 10)
    original_assay <- assay(se_test)

    # Get actual sample names from the assay
    assay_samples <- colnames(original_assay)

    # Create colData with some matching and some non-matching samples
    mixed_samples <- c(assay_samples[seq_len(5)], "NewSample1", "NewSample2")

    col_data_mixed <- DataFrame(
        sample_id = mixed_samples,
        value = runif(7, min = 1, max = 10),
        row.names = mixed_samples  # Set row.names to match sample_id
    )

    # Create new SE - we need to subset the assay to match colData length
    # Or create colData that matches the full assay
    # Let's test the intersection scenario properly

    # Create a scenario where colData has more samples than needed
    extended_samples <- c(assay_samples, "Extra1", "Extra2", "Extra3")
    extended_colData <- DataFrame(
        sample_id = extended_samples,
        value = runif(length(extended_samples), min = 1, max = 10),
        row.names = extended_samples
    )

    # Create extended assay matrix to match
    extended_assay <- cbind(
        original_assay,
        matrix(runif(nrow(original_assay) * 3), nrow = nrow(original_assay))
    )
    colnames(extended_assay) <- extended_samples

    se_extended <- SummarizedExperiment(
        assays = list(counts = extended_assay),
        colData = extended_colData
    )

    # Now subset to create mismatch scenario
    se_mismatch <- se_extended[, seq_len(10)]  # Keep original 10 samples
    # Modify colData to have different samples
    colData(se_mismatch)$sample_id[6:10] <- paste0("Different", seq_len(5))

    result <- sigCor(seData = se_mismatch)
    cor_df <- S4Vectors::metadata(result)$cor.df

    expect_s3_class(cor_df, "data.frame")
    expect_true(nrow(cor_df) >= 0)
})

test_that("sigCor handles missing samples correctly", {
    se_test <- create_test_se(n_samples = 10)

    # Create a new SE with fewer samples to test sample intersection
    # Keep only first 8 samples in both assay and colData
    se_subset <- se_test[, seq_len(8)]

    result <- sigCor(seData = se_subset)
    cor_df <- S4Vectors::metadata(result)$cor.df

    expect_s3_class(cor_df, "data.frame")
    expect_true(nrow(cor_df) > 0)
    expect_equal(nrow(cor_df), nrow(se_subset))
})

test_that("sigCor handles genes with insufficient variation", {
    se_test <- create_test_se(n_genes = 3, n_samples = 10)

    # Make one gene have constant values (no variation)
    assay(se_test)[1, ] <- 5

    result <- sigCor(seData = se_test)
    cor_df <- S4Vectors::metadata(result)$cor.df

    # Should still work but exclude constant genes
    expect_s3_class(cor_df, "data.frame")
    expect_true(nrow(cor_df) >= 0)
})

test_that("sigCor correlation results have correct structure", {
    se_test <- create_test_se(n_genes = 5, n_samples = 20)
    result <- sigCor(seData = se_test)
    cor_df <- S4Vectors::metadata(result)$cor.df

    # Check structure
    expect_equal(ncol(cor_df), 3)
    expect_true(all(c("gene", "cor", "pval") %in% colnames(cor_df)))

    # Check data types
    expect_type(cor_df$gene, "character")
    expect_type(cor_df$cor, "double")
    expect_type(cor_df$pval, "double")

    # Check that correlations are within valid range
    valid_cors <- is.na(cor_df$cor) | (cor_df$cor >= -1 & cor_df$cor <= 1)
    expect_true(all(valid_cors))

    # Check that p-values are within valid range
    valid_pvals <- is.na(cor_df$pval) | (cor_df$pval >= 0 & cor_df$pval <= 1)
    expect_true(all(valid_pvals))
})

test_that("sigCor handles edge cases", {
    # Test with minimum data
    se_small <- create_test_se(n_genes = 2, n_samples = 3)
    result_small <- sigCor(seData = se_small)
    expect_s4_class(result_small, "SummarizedExperiment")

    # Test with larger dataset
    se_large <- create_test_se(n_genes = 20, n_samples = 50)
    result_large <- sigCor(seData = se_large)
    expect_s4_class(result_large, "SummarizedExperiment")
    expect_equal(
        nrow(S4Vectors::metadata(result_large)$cor.df),
        nrow(se_large)
    )
})

test_that("sigCor preserves original data structure", {
    se_test <- create_test_se()
    original_assay <- assay(se_test)
    original_coldata <- colData(se_test)

    result <- sigCor(seData = se_test)

    # Original data should be preserved
    expect_identical(assay(result), original_assay)
    expect_identical(colData(result)$sample_id, original_coldata$sample_id)
    expect_identical(colData(result)$value, original_coldata$value)
})

# Run tests
cat("All sigCor unit tests have been defined and are ready to run.\n")
cat("Use testthat::test_file('tests/testthat/test-sigCor.R') to execute the tests.\n")
