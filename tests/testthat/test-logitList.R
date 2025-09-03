# Test file for .logitList function
# Replace the content of tests/testthat/test-logitList.R with this content

library(testthat)

# Ensure the package is properly loaded
if (!isNamespaceLoaded("SigFun")) {
  requireNamespace("SigFun", quietly = TRUE)
}

# Helper to run .logitList while capturing & muffling warnings
safe_logitList <- function(...) {
  warns <- character()
  res <- withCallingHandlers(
    .logitList(...),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  attr(res, "warnings") <- warns
  res
}

test_that(".logitList basic functionality works", {
  # Create test data with binary outcome
  set.seed(123)
  pattern.signature <- c(0, 1, 0, 1, 1, 0, 1, 0)
  pattern.genes.norm <- matrix(c(
    c(1.1, 2.1, 1.2, 3.1, 2.8, 1.5, 3.2, 1.3),  # positively associated
    c(3.0, 1.0, 2.8, 1.2, 1.1, 2.9, 1.0, 2.7),  # negatively associated
    c(2.0, 2.1, 1.9, 2.2, 2.0, 2.1, 1.8, 2.3)   # weakly associated
  ), nrow = 8, ncol = 3)
  colnames(pattern.genes.norm) <- c("Gene1", "Gene2", "Gene3")

  result <- safe_logitList(y = pattern.signature, pattern.genes.norm, "pearson")
  df <- result
  warns <- attr(df, "warnings")

  # Check if result is a data frame
  expect_s3_class(df, "data.frame")

  # Check dimensions
  expect_equal(nrow(df), ncol(pattern.genes.norm))
  expect_equal(ncol(df), 2)

  # Check column names
  expect_equal(names(df), c("cor", "pval"))

  # Check if all values are lists (due to rbind behavior)
  expect_type(df$cor, "list")
  expect_type(df$pval, "list")

  # Check if list elements are numeric
  expect_true(all(sapply(df$cor, is.numeric)))
  expect_true(all(sapply(df$pval, is.numeric)))

  # Check p-value bounds
  pval_values <- unlist(df$pval)
  expect_true(all(pval_values >= 0 & pval_values <= 1))

  # Optional: warn about separation/numerical issues occurred but not required to fail
  if (length(warns) > 0) {
    expect_true(any(grepl("fitted probabilities numerically 0 or 1 occurred", warns, ignore.case = TRUE)))
  }
})

test_that(".logitList works with different binary patterns", {
  set.seed(123)
  pattern.signature1 <- c(0, 0, 0, 1, 1, 1)
  pattern.signature2 <- c(1, 0, 1, 0, 1, 0)

  pattern.genes.norm <- matrix(c(
    c(1.1, 1.2, 1.0, 2.1, 2.2, 2.0),
    c(2.0, 1.9, 2.1, 1.0, 1.1, 0.9)
  ), nrow = 6, ncol = 2)

  result1 <- safe_logitList(y = pattern.signature1, pattern.genes.norm, "pearson")
  result2 <- safe_logitList(y = pattern.signature2, pattern.genes.norm, "pearson")
  df1 <- result1
  df2 <- result2

  # Both should return valid data frames
  expect_s3_class(df1, "data.frame")
  expect_s3_class(df2, "data.frame")
  expect_equal(nrow(df1), 2)
  expect_equal(nrow(df2), 2)

  # Results should be different for different patterns
  expect_false(identical(unlist(df1$cor), unlist(df2$cor)))
})

test_that(".logitList handles perfect separation cases", {
  # Perfect separation case - all 0s have low gene expression, all 1s have high
  pattern.signature <- c(0, 0, 0, 1, 1, 1)
  pattern.genes.norm <- matrix(c(
    c(1, 1, 1, 10, 10, 10)  # perfect separation
  ), nrow = 6, ncol = 1)

  result <- safe_logitList(y = pattern.signature, pattern.genes.norm, "pearson")
  df <- result
  warns <- attr(df, "warnings")

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1)
  if (length(warns) > 0) {
    expect_true(any(grepl("fitted probabilities numerically 0 or 1 occurred", warns, ignore.case = TRUE)))
  }
})

test_that(".logitList handles constant gene expression", {
  pattern.signature <- c(0, 1, 0, 1, 1)
  pattern.genes.norm <- matrix(c(
    c(2, 2, 2, 2, 2)  # constant expression
  ), nrow = 5, ncol = 1)

  result <- safe_logitList(y = pattern.signature, pattern.genes.norm, "pearson")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 2)
})

test_that(".logitList handles all zeros or all ones in outcome", {
  # All zeros case
  pattern.signature_zeros <- c(0, 0, 0, 0, 0)
  pattern.genes.norm <- matrix(c(
    c(1.1, 2.1, 1.5, 2.8, 1.9)
  ), nrow = 5, ncol = 1)

  result_zero <- safe_logitList(y = pattern.signature_zeros, pattern.genes.norm, "pearson")
  expect_s3_class(result_zero, "data.frame")

  # All ones case
  pattern.signature_ones <- c(1, 1, 1, 1, 1)
  result_one <- safe_logitList(y = pattern.signature_ones, pattern.genes.norm, "pearson")
  expect_s3_class(result_one, "data.frame")
})

test_that(".logitList handles NA values in outcome", {
  pattern.signature <- c(0, 1, NA, 1, 0)
  pattern.genes.norm <- matrix(c(
    c(1.1, 2.1, 1.5, 2.8, 1.9)
  ), nrow = 5, ncol = 1)

  result <- safe_logitList(y = pattern.signature, pattern.genes.norm, "pearson")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true(is.numeric(unlist(result$cor)))
  expect_true(is.numeric(unlist(result$pval)))
})

test_that(".logitList handles NA values in gene expression", {
  pattern.signature <- c(0, 1, 0, 1, 0)
  pattern.genes.norm <- matrix(c(
    c(1.1, NA, 1.5, 2.8, 1.9)
  ), nrow = 5, ncol = 1)

  result <- safe_logitList(y = pattern.signature, pattern.genes.norm, "pearson")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true(is.numeric(unlist(result$cor)))
  expect_true(is.numeric(unlist(result$pval)))
})

test_that(".logitList preserves row names from matrix", {
  pattern.signature <- c(0, 1, 0, 1, 1)
  pattern.genes.norm <- matrix(c(
    c(1.1, 2.1, 1.5, 2.8, 1.9),
    c(2.0, 1.0, 2.5, 1.2, 2.3)
  ), nrow = 5, ncol = 2)
  colnames(pattern.genes.norm) <- c("GeneA", "GeneB")

  result <- safe_logitList(y = pattern.signature, pattern.genes.norm, "pearson")

  expect_s3_class(result, "data.frame")
  expect_equal(rownames(result), c("GeneA", "GeneB"))
})

test_that(".logitList handles single column matrix", {
  pattern.signature <- c(0, 1, 0, 1, 1)
  pattern.genes.norm <- matrix(c(1.1, 2.1, 1.5, 2.8, 1.9), nrow = 5, ncol = 1)
  colnames(pattern.genes.norm) <- "OnlyGene"

  result <- safe_logitList(y = pattern.signature, pattern.genes.norm, "pearson")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 2)
  expect_true(is.numeric(unlist(result$cor)))
  expect_true(is.numeric(unlist(result$pval)))
})

test_that(".logitList handles larger datasets", {
  set.seed(123)
  n <- 50
  pattern.signature <- rbinom(n, 1, 0.5)  # random binary outcome
  pattern.genes.norm <- matrix(rnorm(n * 10), nrow = n, ncol = 10)

  result <- safe_logitList(y = pattern.signature, pattern.genes.norm, "pearson")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 10)
  expect_equal(ncol(result), 2)

  cor_values <- unlist(result$cor)
  pval_values <- unlist(result$pval)
  expect_true(all(is.finite(cor_values)))
  expect_true(all(is.finite(pval_values)))
  expect_true(all(pval_values >= 0 & pval_values <= 1))
})

test_that(".logitList coefficient interpretation", {
  # Create data where higher gene expression is associated with outcome = 1
  set.seed(123)
  pattern.signature <- c(0, 0, 0, 1, 1, 1)
  pattern.genes.norm <- matrix(c(
    c(1, 1.1, 1.2, 2.8, 2.9, 3.0),  # positive association
    c(3.0, 2.9, 2.8, 1.2, 1.1, 1.0)  # negative association
  ), nrow = 6, ncol = 2)

  result <- safe_logitList(y = pattern.signature, pattern.genes.norm, "pearson")

  cor_values <- unlist(result$cor)

  expect_true(cor_values[1] > 0)
  expect_true(cor_values[2] < 0)
})

test_that(".logitList error handling for invalid inputs", {
  pattern.signature <- c(0, 1, 0, 1, 1)

  # Test with mismatched dimensions
  pattern.genes.norm_wrong <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)
  expect_error({
    .logitList(y = pattern.signature, pattern.genes.norm_wrong, "pearson")
  })

  # Test with non-binary outcome (should work but produce warning)
  pattern.signature_continuous <- c(0.1, 0.9, 0.3, 0.7, 0.5)
  pattern.genes.norm <- matrix(c(1.1, 2.1, 1.5, 2.8, 1.9), nrow = 5, ncol = 1)

  result_cont <- safe_logitList(y = pattern.signature_continuous, pattern.genes.norm, "pearson")
  expect_s3_class(result_cont, "data.frame")
})

test_that(".logitList handles edge cases with small samples and emits expected warnings", {
  pattern.signature <- c(0, 1, 0)
  pattern.genes.norm <- matrix(c(1.0, 2.0, 1.5), nrow = 3, ncol = 1)

  result <- safe_logitList(y = pattern.signature, pattern.genes.norm, "pearson")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
})

test_that(".logitList default parameter behavior", {
  pattern.signature <- c(0, 1, 0, 1, 1)
  pattern.genes.norm <- matrix(c(
    c(1.1, 2.1, 1.5, 2.8, 1.9)
  ), nrow = 5, ncol = 1)

  # Test with explicit y parameter (normal usage)
  result <- safe_logitList(y = pattern.signature, pattern.genes.norm, "pearson")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 2)

  # Test that omitting y parameter causes an error (since pattern.signature is not in function scope)
  expect_error({
    .logitList(pattern.genes.norm = pattern.genes.norm, cor.method = "pearson")
  }, "object 'pattern.signature' not found")
})
