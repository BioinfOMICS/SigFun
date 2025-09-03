# Test file for .checkSE function
# Replace the content of tests/testthat/test-checkSE.R with this content

library(testthat)
library(SummarizedExperiment)

# Ensure the package namespace is loaded
if (!isNamespaceLoaded("SigFun")) {
  requireNamespace("SigFun", quietly = TRUE)
}

test_that(".checkSE succeeds for a valid SummarizedExperiment", {
  se <- SummarizedExperiment(
    assays = list(dummy = matrix(1:4, nrow = 2)),
    metadata = list(gseaResult = data.frame())
  )
  expect_silent(
    .checkSE(se)
  )
})

test_that(".checkSE errors when object is not a SummarizedExperiment", {
  expect_error(
    .checkSE(42),
    regexp = "inherits\\(se,.*SummarizedExperiment.*\\) is not TRUE"
  )
})

test_that(".checkSE errors when metadata slot exists but no gseaResult", {
  se <- SummarizedExperiment(
    assays = list(dummy = matrix(1:4, nrow = 2)),
    metadata = list(other = "foo")
  )
  expect_error(
    .checkSE(se),
    regexp = "gseaResult.*is not TRUE"
  )
})
