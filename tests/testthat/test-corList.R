# Test file for .corList function
# Replace the content of tests/testthat/test-corList.R with this content

library(testthat)

# Ensure the package is properly loaded
if (!isNamespaceLoaded("SigFun")) {
  requireNamespace("SigFun", quietly = TRUE)
}

test_that(".corList basic functionality works", {
  # Create test data
  set.seed(123)
  pattern.signature <- c(1, 2, 3, 4, 5)
  pattern.genes.norm <- matrix(c(
    c(1.1, 2.1, 2.9, 4.1, 5.0),  # positively correlated
    c(5.0, 4.0, 3.0, 2.0, 1.0),  # negatively correlated
    c(2.0, 3.0, 1.0, 5.0, 4.0)   # mixed correlation
  ), nrow = 5, ncol = 3)
  colnames(pattern.genes.norm) <- c("Gene1", "Gene2", "Gene3")

  result <- .corList(pattern.signature, pattern.genes.norm, "pearson")

  # Check if result is a data frame
  expect_s3_class(result, "data.frame")

  # Check dimensions
  expect_equal(nrow(result), ncol(pattern.genes.norm))
  expect_equal(ncol(result), 2)

  # Check column names
  expect_equal(names(result), c("cor", "pval"))

  # Check if all values are numeric (they come as lists from rbind)
  expect_type(result$cor, "list")
  expect_type(result$pval, "list")

  # Check if list elements are numeric
  expect_true(all(sapply(result$cor, is.numeric)))
  expect_true(all(sapply(result$pval, is.numeric)))

  # Check correlation bounds (extract numeric values from lists)
  cor_values <- unlist(result$cor)
  pval_values <- unlist(result$pval)
  expect_true(all(cor_values >= -1 & cor_values <= 1))

  # Check p-value bounds
  expect_true(all(pval_values >= 0 & pval_values <= 1))
})

test_that(".corList works with different correlation methods", {
  set.seed(123)
  pattern.signature <- c(1, 2, 3, 4, 5)
  pattern.genes.norm <- matrix(c(
    c(1.1, 2.1, 2.9, 4.1, 5.0),
    c(5.0, 4.0, 3.0, 2.0, 1.0)
  ), nrow = 5, ncol = 2)

  # Test Pearson correlation
  result_pearson <- SigFun:::.corList(pattern.signature, pattern.genes.norm, "pearson")
  expect_s3_class(result_pearson, "data.frame")
  expect_equal(nrow(result_pearson), 2)

  # Test Spearman correlation
  result_spearman <- SigFun:::.corList(pattern.signature, pattern.genes.norm, "spearman")
  expect_s3_class(result_spearman, "data.frame")
  expect_equal(nrow(result_spearman), 2)

  # Test Kendall correlation
  result_kendall <- SigFun:::.corList(pattern.signature, pattern.genes.norm, "kendall")
  expect_s3_class(result_kendall, "data.frame")
  expect_equal(nrow(result_kendall), 2)

  # Results should be different for different methods
  expect_false(identical(unlist(result_pearson$cor), unlist(result_spearman$cor)))
})

test_that(".corList handles perfect correlations", {
  # Perfect positive correlation
  pattern.signature <- c(1, 2, 3, 4, 5)
  pattern.genes.norm <- matrix(c(
    c(2, 4, 6, 8, 10),  # perfect positive correlation (r = 1)
    c(5, 4, 3, 2, 1)    # perfect negative correlation (r = -1)
  ), nrow = 5, ncol = 2)

  result <- SigFun:::.corList(pattern.signature, pattern.genes.norm, "pearson")

  # Check perfect correlations (extract from lists)
  expect_equal(unlist(result$cor)[1], 1, tolerance = 1e-10)
  expect_equal(unlist(result$cor)[2], -1, tolerance = 1e-10)

  # Check p-values for perfect correlations (should be very small)
  expect_true(unlist(result$pval)[1] < 0.01)
  expect_true(unlist(result$pval)[2] < 0.01)
})

test_that(".corList handles constant values and captures expected warning without polluting summary", {
  # Constant signature and one constant gene expression
  pattern.signature <- c(3, 3, 3, 3, 3)
  pattern.genes.norm <- matrix(c(
    c(1, 2, 3, 4, 5),
    c(2, 2, 2, 2, 2)  # constant gene expression
  ), nrow = 5, ncol = 2)

  warns <- character()
  result <- withCallingHandlers(
    SigFun:::.corList(pattern.signature, pattern.genes.norm, "pearson"),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")  # suppress it from bubbling up
    }
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)

  # Confirm the expected warning about zero standard deviation occurred
  expect_true(any(grepl("standard deviation is zero", warns, ignore.case = TRUE)))
})



test_that(".corList handles NA values correctly", {
  pattern.signature <- c(1, 2, NA, 4, 5)
  pattern.genes.norm <- matrix(c(
    c(1.1, 2.1, 2.9, 4.1, 5.0),
    c(5.0, NA, 3.0, 2.0, 1.0)
  ), nrow = 5, ncol = 2)

  result <- SigFun:::.corList(pattern.signature, pattern.genes.norm, "pearson")

  # Should handle pairwise complete observations
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)

  # Results should not be NA (due to pairwise.complete.obs)
  expect_false(is.na(unlist(result$cor)[1]))
  expect_false(is.na(unlist(result$pval)[1]))
})

test_that(".corList handles single column matrix", {
  pattern.signature <- c(1, 2, 3, 4, 5)
  pattern.genes.norm <- matrix(c(1.1, 2.1, 2.9, 4.1, 5.0), nrow = 5, ncol = 1)
  colnames(pattern.genes.norm) <- "OnlyGene"

  result <- SigFun:::.corList(pattern.signature, pattern.genes.norm, "pearson")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 2)
  expect_true(abs(unlist(result$cor)[1]) <= 1)
})

test_that(".corList preserves row names from matrix", {
  pattern.signature <- c(1, 2, 3, 4, 5)
  pattern.genes.norm <- matrix(c(
    c(1.1, 2.1, 2.9, 4.1, 5.0),
    c(5.0, 4.0, 3.0, 2.0, 1.0)
  ), nrow = 5, ncol = 2)
  colnames(pattern.genes.norm) <- c("GeneA", "GeneB")

  result <- SigFun:::.corList(pattern.signature, pattern.genes.norm, "pearson")

  # Check if row names are preserved from column names
  expect_equal(rownames(result), c("GeneA", "GeneB"))
})

test_that(".corList handles large datasets", {
  set.seed(123)
  pattern.signature <- rnorm(100)
  pattern.genes.norm <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)

  result <- SigFun:::.corList(pattern.signature, pattern.genes.norm, "pearson")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 50)
  expect_equal(ncol(result), 2)
  expect_true(all(unlist(result$cor) >= -1 & unlist(result$cor) <= 1, na.rm = TRUE))
  expect_true(all(unlist(result$pval) >= 0 & unlist(result$pval) <= 1, na.rm = TRUE))
})

test_that(".corList handles edge cases with very small samples", {
  # Test with minimum sample size (n=3)
  pattern.signature <- c(1, 2, 3)
  pattern.genes.norm <- matrix(c(
    c(1.1, 2.1, 2.9),
    c(3.0, 2.0, 1.0)
  ), nrow = 3, ncol = 2)

  result <- SigFun:::.corList(pattern.signature, pattern.genes.norm, "pearson")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(unlist(result$cor) >= -1 & unlist(result$cor) <= 1, na.rm = TRUE))
})

test_that(".corList error handling for invalid inputs", {
  pattern.signature <- c(1, 2, 3, 4, 5)

  # Test with mismatched dimensions
  pattern.genes.norm_wrong <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)
  expect_error(SigFun:::.corList(pattern.signature, pattern.genes.norm_wrong, "pearson"))

  # Test with invalid correlation method
  pattern.genes.norm <- matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1)
  expect_error(SigFun:::.corList(pattern.signature, pattern.genes.norm, "invalid_method"))
})

test_that(".corList handles all NA values in a column", {
  pattern.signature <- c(1, 2, 3, 4, 5)
  pattern.genes.norm <- matrix(c(
    c(1.1, 2.1, 2.9, 4.1, 5.0),
    c(NA, NA, NA, NA, NA)  # all NA column
  ), nrow = 5, ncol = 2)

  # This should produce an error for the all-NA column
  expect_error({
    result <- SigFun:::.corList(pattern.signature, pattern.genes.norm, "pearson")
  }, "not enough finite observations")
})

test_that(".corList numerical stability", {
  # Test with very similar but not identical values
  pattern.signature <- c(1.0000001, 1.0000002, 1.0000003, 1.0000004, 1.0000005)
  pattern.genes.norm <- matrix(c(
    c(2.0000001, 2.0000002, 2.0000003, 2.0000004, 2.0000005)
  ), nrow = 5, ncol = 1)

  result <- SigFun:::.corList(pattern.signature, pattern.genes.norm, "pearson")

  expect_s3_class(result, "data.frame")
  expect_true(abs(unlist(result$cor)[1] - 1) < 1e-10)  # Should be very close to 1
  expect_true(unlist(result$pval)[1] < 0.05)  # Should be significant
})
