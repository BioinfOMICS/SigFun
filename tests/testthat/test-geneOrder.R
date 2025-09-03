# Test file for .geneOrder function
# Replace the content of tests/testthat/test-geneOrder.R with this content

library(testthat)
# Ensure the SigFun namespace is loaded without attaching
if (!isNamespaceLoaded("SigFun")) requireNamespace("SigFun", quietly = TRUE)

test_that(".geneOrder returns a character vector of unique genes", {
  df <- data.frame(
    categoryID = rep(c("A", "B"), each = 3),
    Gene       = rep(c("g1", "g2", "g3"), times = 2),
    Coef       = c(0, 1, 2, 0, 2, 1),
    stringsAsFactors = FALSE
  )
  out <- .geneOrder(df, distfun = "euclidean", hclustfun = "complete")
  expect_type(out, "character")
  expect_equal(length(out), 3)
  expect_setequal(out, c("g1", "g2", "g3"))
})

test_that(".geneOrder is deterministic with NA values replaced by zero", {
  df_na <- data.frame(
    categoryID = rep(c("A", "B"), each = 3),
    Gene       = rep(c("g1", "g2", "g3"), times = 2),
    Coef       = c(0, 1, 2, NA, 2, 1),
    stringsAsFactors = FALSE
  )
  out_na <- .geneOrder(df_na, distfun = "euclidean", hclustfun = "complete")
  df <- data.frame(
    categoryID = rep(c("A", "B"), each = 3),
    Gene       = rep(c("g1", "g2", "g3"), times = 2),
    Coef       = c(0, 1, 2, 0, 2, 1),
    stringsAsFactors = FALSE
  )
  out <- .geneOrder(df, distfun = "euclidean", hclustfun = "complete")
  expect_equal(out_na, out)
})

test_that(".geneOrder errors when fewer than two genes are present", {
  df_single <- data.frame(
    categoryID = "A",
    Gene       = "g1",
    Coef       = 5,
    stringsAsFactors = FALSE
  )
  expect_error(
    .geneOrder(df_single, distfun = "euclidean", hclustfun = "complete"),
    regexp = "must have n >= 2"
  )
})
