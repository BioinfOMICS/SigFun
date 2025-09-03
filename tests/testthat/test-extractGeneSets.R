# Test file for .extractGeneSets function
# Replace the content of tests/testthat/test-extractGeneSets.R with this content

library(testthat)

if (!isNamespaceLoaded("SigFun"))
  requireNamespace("SigFun", quietly = TRUE)

# Define a minimal dummy S4 class to simulate gseaResult
setClass("dummyGsea",
         slots = c(result = "data.frame"),
         package = "SigFun")

make_dummyGsea <- function(ids, descs) {
  new("dummyGsea",
      result = data.frame(
        ID          = ids,
        Description = descs,
        stringsAsFactors = FALSE
      ))
}

test_that(".extractGeneSets subsets correctly for numeric showCategory", {
  gr       <- make_dummyGsea(ids = c("A", "B", "C"), descs = c("D1", "D2", "D3"))
  geneSets <- list(A = c("gA1", "gA2"), B = c("gB"), C = c("gC"))
  expected <- geneSets[1:2]
  names(expected) <- gr@result$Description[1:2]

  # Expected output should match the actual structure returned by .list2df
  expected_df <- data.frame(
    categoryID = c("D1", "D1", "D2"),
    Gene = c("gA1", "gA2", "gB"),
    stringsAsFactors = FALSE
  )

  # Stub DOSE::geneInCategory
  orig_gic <- get("geneInCategory", envir = asNamespace("DOSE"))
  unlockBinding("geneInCategory", asNamespace("DOSE"))
  assign("geneInCategory", function(x) {
    expect_identical(x, gr)
    geneSets
  }, envir = asNamespace("DOSE"))
  lockBinding("geneInCategory", asNamespace("DOSE"))

  # Don't stub .list2df - let it work naturally
  on.exit({
    unlockBinding("geneInCategory", asNamespace("DOSE"))
    assign("geneInCategory", orig_gic, envir = asNamespace("DOSE"))
    lockBinding("geneInCategory", asNamespace("DOSE"))
  }, add = TRUE)

  out <- .extractGeneSets(gr, showCategory = 2)
  expect_identical(out, expected_df)
})

test_that(".extractGeneSets subsets correctly for character showCategory", {
  gr       <- make_dummyGsea(ids = c("X", "Y", "Z"), descs = c("DescX", "DescY", "DescZ"))
  geneSets <- list(X = c("gX1"), Y = c("gY1", "gY2"), Z = c("gZ1"))
  showChar <- c("DescZ", "DescX")

  # Expected output based on the selected categories
  expected_df <- data.frame(
    categoryID = c("DescZ", "DescX"),
    Gene = c("gZ1", "gX1"),
    stringsAsFactors = FALSE
  )

  # Stub DOSE::geneInCategory
  orig_gic <- get("geneInCategory", envir = asNamespace("DOSE"))
  unlockBinding("geneInCategory", asNamespace("DOSE"))
  assign("geneInCategory", function(x) {
    expect_identical(x, gr)
    geneSets
  }, envir = asNamespace("DOSE"))
  lockBinding("geneInCategory", asNamespace("DOSE"))

  on.exit({
    unlockBinding("geneInCategory", asNamespace("DOSE"))
    assign("geneInCategory", orig_gic, envir = asNamespace("DOSE"))
    lockBinding("geneInCategory", asNamespace("DOSE"))
  }, add = TRUE)

  out <- .extractGeneSets(gr, showCategory = showChar)
  expect_identical(out, expected_df)
})

test_that(".extractGeneSets errors when character showCategory not found", {
  gr       <- make_dummyGsea(ids = c("M"), descs = c("DescM"))
  geneSets <- list(M = c("gm"))

  # Stub DOSE::geneInCategory
  orig_gic <- get("geneInCategory", envir = asNamespace("DOSE"))
  unlockBinding("geneInCategory", asNamespace("DOSE"))
  assign("geneInCategory", function(x) geneSets, envir = asNamespace("DOSE"))
  lockBinding("geneInCategory", asNamespace("DOSE"))

  on.exit({
    unlockBinding("geneInCategory", asNamespace("DOSE"))
    assign("geneInCategory", orig_gic, envir = asNamespace("DOSE"))
    lockBinding("geneInCategory", asNamespace("DOSE"))
  }, add = TRUE)

  expect_error(
    .extractGeneSets(gr, showCategory = "NotThere"),
    regexp = "cannot be found"
  )
})
