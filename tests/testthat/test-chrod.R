# Test file for .chrod function
# Replace the content of tests/testthat/test-chrod.R with this content

library(testthat)
library(circlize)

# Ensure the SigFun namespace is loaded without attaching
if (!isNamespaceLoaded("SigFun")) requireNamespace("SigFun", quietly = TRUE)

test_that(".chrod runs without error and clears circos state", {
  # Prepare a minimal valid data frame: from, to, value
  df <- data.frame(
    name  = c("P1", "P2", "P1"),
    Gene  = c("G1", "G2", "G3"),
    value = c(1, 1, 1),
    stringsAsFactors = FALSE
  )
  fontSize <- 1.2

  # Clear any existing circos state
  circlize::circos.clear()

  # Should not error and should return NULL invisibly
  expect_silent({
    res <- .chrod(df, fontSize)
    expect_null(res)
  })

  # After running, circos should have been cleared again
  info <- circlize::circos.info(plot = FALSE)
  # No sectors remain
  expect_true(length(info$sector.index) == 0)
})
