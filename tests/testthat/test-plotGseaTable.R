# Test file for .plotGseaTable function
# Replace the content of tests/testthat/test-plotGseaTable.R with this content

library(testthat)
library(ggplot2)
library(dplyr)
library(cowplot)

# Ensure the package namespace is loaded
if (!isNamespaceLoaded("SigFun")) {
  requireNamespace("SigFun", quietly = TRUE)
}

test_that(".plotGseaTable basic functionality returns ggplot and handles deprecation warning", {
  # Prepare test data: one valid pathway and one invalid (should be filtered out)
  stats <- c(geneA = 2, geneB = -1, geneC = 1.5, geneD = -0.5)
  Pathways <- list(
    valid_path = c("geneA", "geneC"),
    bad_path   = c("nonexistent")
  )
  fgseaRes <- data.frame(
    ID     = c("valid_path", "bad_path"),
    NES    = c(1.5, -2.0),
    pvalue = c(0.01, 0.05),
    stringsAsFactors = FALSE
  )

  # Expect a deprecation warning when using the render argument
  expect_warning(
    p <- .plotGseaTable(
      pathways              = Pathways,
      stats                 = stats,
      fgseaRes              = fgseaRes,
      gseaParam             = 1,
      pathwayLabelStyle     = list(),
      render                = TRUE,  # triggers deprecation warning
      headerLabelStyle      = list(),
      valueStyle            = list(),
      axisLabelStyle        = list()
    ),
    regexp = "render argument is deprecated"
  )

  # Should return a non-null ggplot object
  expect_s3_class(p, "ggplot")
  expect_false(is.null(p))
})

test_that(".plotGseaTable handles case with all pathways valid", {
  # Prepare data where all pathways are valid
  stats <- c(g1 = 0.5, g2 = -0.2, g3 = 1.2, g4 = -0.8)
  Pathways <- list(
    p_alpha = c("g1", "g3"),
    p_beta  = c("g2", "g4")
  )
  fgseaRes <- data.frame(
    ID     = c("p_alpha", "p_beta"),
    NES    = c(0.8, -1.1),
    pvalue = c(0.02, 0.03),
    stringsAsFactors = FALSE
  )

  # Should not emit any warnings when render is not provided
  expect_no_warning(
    p <- .plotGseaTable(
      pathways              = Pathways,
      stats                 = stats,
      fgseaRes              = fgseaRes,
      gseaParam             = 2,
      pathwayLabelStyle     = list(size = 5),
      headerLabelStyle      = list(),
      valueStyle            = list(),
      axisLabelStyle        = list()
    )
  )

  expect_s3_class(p, "ggplot")
  expect_false(is.null(p))
})

test_that(".plotGseaTable handles case with no valid pathways", {
  # Prepare data where none of the pathways overlap the stats
  stats <- c(geneA = 2, geneB = -1)
  Pathways <- list(
    bad_path1 = c("nonexistent1"),
    bad_path2 = c("nonexistent2")
  )
  fgseaRes <- data.frame(
    ID     = c("bad_path1", "bad_path2"),
    NES    = c(1.5, -2.0),
    pvalue = c(0.01, 0.05),
    stringsAsFactors = FALSE
  )

  # Expect a warning from max() when there are no non-missing arguments
  expect_warning(
    p <- .plotGseaTable(
      pathways              = Pathways,
      stats                 = stats,
      fgseaRes              = fgseaRes,
      gseaParam             = 1,
      pathwayLabelStyle     = list(),
      headerLabelStyle      = list(),
      valueStyle            = list(),
      axisLabelStyle        = list()
    ),
    regexp = "no non-missing arguments to max"
  )

  expect_s3_class(p, "ggplot")
  expect_false(is.null(p))
})

test_that(".plotGseaTable handles mixed valid and invalid pathways", {
  # Prepare data with a mix of valid and invalid pathways
  stats <- c(geneA = 2, geneB = -1, geneC = 1.5)
  Pathways <- list(
    valid_path    = c("geneA", "geneC"),  # valid
    bad_path      = c("nonexistent"),     # invalid
    another_valid = c("geneB")            # valid
  )
  fgseaRes <- data.frame(
    ID     = c("valid_path", "bad_path", "another_valid"),
    NES    = c(1.5, -2.0, 0.8),
    pvalue = c(0.01, 0.05, 0.03),
    stringsAsFactors = FALSE
  )

  expect_no_warning(
    p <- .plotGseaTable(
      pathways              = Pathways,
      stats                 = stats,
      fgseaRes              = fgseaRes,
      gseaParam             = 1,
      pathwayLabelStyle     = list(),
      headerLabelStyle      = list(),
      valueStyle            = list(),
      axisLabelStyle        = list()
    )
  )

  expect_s3_class(p, "ggplot")
  expect_false(is.null(p))
})

test_that(".plotGseaTable parameter validation and edge cases", {
  # Test various gseaParam values and ensure no errors
  stats <- c(geneA = 2, geneB = -1)
  Pathways <- list(valid_path = c("geneA"))
  fgseaRes <- data.frame(
    ID     = "valid_path",
    NES    = 1.5,
    pvalue = 0.01,
    stringsAsFactors = FALSE
  )

  expect_no_error({
    p1 <- .plotGseaTable(
      pathways              = Pathways,
      stats                 = stats,
      fgseaRes              = fgseaRes,
      gseaParam             = 0.5,  # less than 1
      pathwayLabelStyle     = list(),
      headerLabelStyle      = list(),
      valueStyle            = list(),
      axisLabelStyle        = list()
    )
  })

  expect_no_error({
    p2 <- .plotGseaTable(
      pathways              = Pathways,
      stats                 = stats,
      fgseaRes              = fgseaRes,
      gseaParam             = 3,    # greater than 1
      pathwayLabelStyle     = list(),
      headerLabelStyle      = list(),
      valueStyle            = list(),
      axisLabelStyle        = list()
    )
  })

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that(".plotGseaTable style parameters work correctly", {
  # Test that custom style lists do not cause errors
  stats <- c(geneA = 2, geneB = -1, geneC = 0.5)
  Pathways <- list(valid_path = c("geneA"), another_path = c("geneB"))
  fgseaRes <- data.frame(
    ID     = c("valid_path", "another_path"),
    NES    = c(1.5, -1.2),
    pvalue = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  expect_no_error({
    p <- .plotGseaTable(
      pathways              = Pathways,
      stats                 = stats,
      fgseaRes              = fgseaRes,
      gseaParam             = 1,
      pathwayLabelStyle     = list(size = 12, hjust = 0.5),
      headerLabelStyle      = list(size = 14),
      valueStyle            = list(size = 10),
      axisLabelStyle        = list(size = 8)
    )
  })

  expect_s3_class(p, "ggplot")
})

test_that(".plotGseaTable handles custom colwidths", {
  # Test custom column widths argument
  stats <- c(geneA = 2, geneB = -1)
  Pathways <- list(valid_path = c("geneA"))
  fgseaRes <- data.frame(
    ID     = "valid_path",
    NES    = 1.5,
    pvalue = 0.01,
    stringsAsFactors = FALSE
  )

  expect_no_error({
    p <- .plotGseaTable(
      pathways              = Pathways,
      stats                 = stats,
      fgseaRes              = fgseaRes,
      gseaParam             = 1,
      colwidths             = c(1, 1, 1.5, 3, 8),  # custom widths
      pathwayLabelStyle     = list(),
      headerLabelStyle      = list(),
      valueStyle            = list(),
      axisLabelStyle        = list()
    )
  })

  expect_s3_class(p, "ggplot")
})

test_that(".plotGseaTable handles empty stats", {
  # Test behavior when stats vector is empty
  stats <- numeric(0)
  Pathways <- list(valid_path = c("geneA"))
  fgseaRes <- data.frame(
    ID     = "valid_path",
    NES    = 1.5,
    pvalue = 0.01,
    stringsAsFactors = FALSE
  )

  result <- suppressWarnings({
    warning_caught <- FALSE
    tryCatch({
      withCallingHandlers({
        p <- .plotGseaTable(
          pathways              = Pathways,
          stats                 = stats,
          fgseaRes              = fgseaRes,
          gseaParam             = 1,
          pathwayLabelStyle     = list(),
          headerLabelStyle      = list(),
          valueStyle            = list(),
          axisLabelStyle        = list()
        )
      }, warning = function(w) {
        if (grepl("no non-missing arguments to max", w$message)) {
          warning_caught <<- TRUE
        }
      })
      list(plot = p, warning_caught = warning_caught)
    }, error = function(e) {
      list(plot = NULL, warning_caught = warning_caught, error = e)
    })
  })

  # Confirm that the expected warning was caught
  expect_true(result$warning_caught, "Expected warning about 'no non-missing arguments to max'")

  # The function should still return a ggplot
  expect_s3_class(result$plot, "ggplot")
})

test_that(".plotGseaTable handles single pathway", {
  # Test with exactly one pathway
  stats <- c(geneA = 2, geneB = -1, geneC = 0.5)
  Pathways <- list(single_path = c("geneA", "geneB"))
  fgseaRes <- data.frame(
    ID     = "single_path",
    NES    = 1.5,
    pvalue = 0.01,
    stringsAsFactors = FALSE
  )

  expect_no_error({
    p <- .plotGseaTable(
      pathways              = Pathways,
      stats                 = stats,
      fgseaRes              = fgseaRes,
      gseaParam             = 1,
      pathwayLabelStyle     = list(),
      headerLabelStyle      = list(),
      valueStyle            = list(),
      axisLabelStyle        = list()
    )
  })

  expect_s3_class(p, "ggplot")
})

test_that(".plotGseaTable handles NULL style parameters", {
  # Test that passing NULL for style lists uses defaults
  stats <- c(geneA = 2, geneB = -1)
  Pathways <- list(valid_path = c("geneA"))
  fgseaRes <- data.frame(
    ID     = "valid_path",
    NES    = 1.5,
    pvalue = 0.01,
    stringsAsFactors = FALSE
  )

  expect_no_error({
    p <- .plotGseaTable(
      pathways              = Pathways,
      stats                 = stats,
      fgseaRes              = fgseaRes,
      gseaParam             = 1,
      pathwayLabelStyle     = NULL,
      headerLabelStyle      = NULL,
      valueStyle            = NULL,
      axisLabelStyle        = NULL
    )
  })

  expect_s3_class(p, "ggplot")
})
