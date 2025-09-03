# Test file for .plotEnrichment function
# Replace the content of tests/testthat/test-plotEnrichment.R with this content

library(testthat)
library(ggplot2)

# Ensure the package namespace is loaded
if (!isNamespaceLoaded("SigFun")) {
  requireNamespace("SigFun", quietly = TRUE)
}

# More robust mock helper: reliably override and restore .plotEnrichmentData
with_mocked_enrichment_data <- function(code) {
  ns <- asNamespace("SigFun")
  has_orig <- exists(".plotEnrichmentData", envir = ns, inherits = FALSE)
  orig <- if (has_orig) get(".plotEnrichmentData", envir = ns) else NULL

  mock <- function(pathway, stats, gseaParam = 1) {
    curve <- data.frame(rank = 1:4, dummy = c(1, 2, 3, 4))
    ticks <- data.frame(rank = 1:4, spreadES = c(0.2, 0.4, 0.6, 0.8))
    list(curve = curve, ticks = ticks)
  }

  if (bindingIsLocked(".plotEnrichmentData", ns)) unlockBinding(".plotEnrichmentData", ns)
  assign(".plotEnrichmentData", mock, envir = ns)
  lockBinding(".plotEnrichmentData", ns)

  on.exit({
    if (bindingIsLocked(".plotEnrichmentData", ns)) unlockBinding(".plotEnrichmentData", ns)
    if (has_orig) {
      assign(".plotEnrichmentData", orig, envir = ns)
    } else {
      if (exists(".plotEnrichmentData", envir = ns, inherits = FALSE)) {
        rm(".plotEnrichmentData", envir = ns)
      }
    }
    lockBinding(".plotEnrichmentData", ns)
  }, add = TRUE)

  force(code)
}

test_that(".plotEnrichment basic plot builds without error and returns ggplot (subset pathway)", {
  with_mocked_enrichment_data({
    # pathway is a subset of stats to avoid full-gene underlying error
    p <- .plotEnrichment(
      pathway = c("gene1"),
      stats = c(gene1 = 1, gene2 = 2),
      color = 1
    )
    expect_s3_class(p, "ggplot")
  })
})

test_that(".plotEnrichment positive color produces '#FD7C69' for tick segments", {
  with_mocked_enrichment_data({
    # pathway is still a subset (avoiding full-gene)
    p <- .plotEnrichment(
      pathway = c("g"),
      stats = c(g = 1, h = 0),
      color = 5,       # positive
      gseaParam = 2,
      ticksSize = 0.3
    )
    built <- ggplot2::ggplot_build(p)
    seg_data <- built$data[[1]]
    expect_true(all(tolower(seg_data$colour) == "#fd7c69"))
    expect_equal(unique(seg_data$linewidth), 0.3)
  })
})

test_that(".plotEnrichment negative color produces '#5F90BB' for tick segments", {
  with_mocked_enrichment_data({
    p <- .plotEnrichment(
      pathway = c("g"),
      stats = c(g = 1, h = 0),
      color = -2  # negative
    )
    built <- ggplot2::ggplot_build(p)
    seg_data <- built$data[[1]]
    expect_true(all(tolower(seg_data$colour) == "#5f90bb"))
  })
})

test_that(".plotEnrichment different ticksSize affects segment linewidth", {
  with_mocked_enrichment_data({
    p_small <- .plotEnrichment(
      pathway = c("g"),
      stats = c(g = 1, h = 0),
      color = 1,
      ticksSize = 0.1
    )
    p_large <- .plotEnrichment(
      pathway = c("g"),
      stats = c(g = 1, h = 0),
      color = 1,
      ticksSize = 0.8
    )
    built_small <- ggplot2::ggplot_build(p_small)
    built_large <- ggplot2::ggplot_build(p_large)
    lw_small <- unique(built_small$data[[1]]$linewidth)
    lw_large <- unique(built_large$data[[1]]$linewidth)
    expect_true(lw_small < lw_large)
  })
})

test_that(".plotEnrichment errors when all genes are selected (no mock)", {
  # Directly trigger the full-gene error case without mocking
  stats <- c(a = 1, b = 2)
  expect_error(
    .plotEnrichment(
      pathway = names(stats),
      stats = stats,
      color = 1
    ),
    regexp = "GSEA statistic is not defined when all genes are selected"
  )
})
