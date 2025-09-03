# Test file for .plotEnrichmentData function
# Replace the content of tests/testthat/test-plotEnrichmentData.R with this content

library(testthat)

# skip if fgsea is not installed, since the function depends on it
skip_if_not_installed("fgsea")

# Ensure the namespace is loaded if it's in SigFun
if (!isNamespaceLoaded("SigFun")) {
  requireNamespace("SigFun", quietly = TRUE)
}

test_that(".plotEnrichmentData basic functionality returns expected structure", {
  stats <- c(a = 3, b = 1, c = -2, d = 0.5)
  pathway <- c("a", "c")
  res <- .plotEnrichmentData(pathway = pathway, stats = stats, gseaParam = 1)

  expected_names <- c("curve", "ticks", "stats", "posES", "negES", "spreadES", "maxAbsStat")
  expect_true(all(expected_names %in% names(res)))

  expect_true(is.data.frame(res$curve))
  expect_true(is.data.frame(res$ticks))
  expect_true(is.data.frame(res$stats))
  expect_true(is.numeric(res$posES))
  expect_true(is.numeric(res$negES))
  expect_true(is.numeric(res$spreadES))
  expect_true(is.numeric(res$maxAbsStat))
  expect_equal(res$spreadES, res$posES - res$negES)

  # Recompute adjusted stats in same way as function does
  ord <- order(rank(-stats))
  statsAdj <- sign(stats[ord]) * (abs(stats[ord]) ^ 1)
  expect_equal(res$maxAbsStat, max(abs(statsAdj)))

  # pathway matching after reordering
  matched <- sort(unique(na.omit(match(pathway, names(statsAdj)))))
  expect_equal(sort(res$ticks$rank), matched)
})

test_that(".plotEnrichmentData respects gseaParam transformation", {
  stats <- c(a = 2, b = -3, c = 1)  # subset pathway so not all genes selected
  pathway <- c("a", "b")
  res2 <- .plotEnrichmentData(pathway = pathway, stats = stats, gseaParam = 2)

  # Extract adjusted stats from returned stats table
  returned_stats <- res2$stats$stat

  # Manually compute expected adjusted stats (with ordering)
  ord <- order(rank(-stats))
  expected_adj <- sign(stats[ord]) * (abs(stats[ord]) ^ 2)

  # Compare as sets, stripping names to avoid name mismatch
  expect_equal(sort(returned_stats), unname(sort(expected_adj)))
})


test_that(".plotEnrichmentData errors when stats contains non-finite values", {
  stats_bad <- c(a = 1, b = Inf, c = 2)
  pathway <- c("a", "c")
  expect_error(
    .plotEnrichmentData(pathway = pathway, stats = stats_bad),
    regexp = "Not all stats values are finite numbers"
  )
})

test_that(".plotEnrichmentData handles pathway with no overlap and captures expected warnings without noisy print", {
  stats <- c(a = 1, b = 2, c = 3)
  pathway <- c("x", "y")  # no overlap

  warns <- character()
  res <- withCallingHandlers(
    .plotEnrichmentData(pathway = pathway, stats = stats),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")  # suppress propagation
    }
  )

  # Expect both kinds of warnings were emitted
  expect_true(any(grepl("no non-missing arguments to max", warns)))
  expect_true(any(grepl("no non-missing arguments to min", warns)))

  # Core structure still present
  expect_true(is.data.frame(res$curve))
  expect_true(is.data.frame(res$ticks))
  expect_equal(length(res$ticks$rank), 0)

  # Current behavior: infinite ES values
  expect_equal(res$posES, -Inf)
  expect_equal(res$negES, Inf)
  expect_true(is.infinite(res$spreadES))
  expect_equal(res$spreadES, res$posES - res$negES)
})
