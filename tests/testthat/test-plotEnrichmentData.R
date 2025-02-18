library(testthat)
library(ggplot2)
library(fgsea)

test_that(".plotEnrichment function works correctly", {
  # Create mock data for testing
  set.seed(123)

  # Simulate a realistic stats vector (sorted)
  mock_stats <- sort(rnorm(200), decreasing = TRUE)
  names(mock_stats) <- paste0("gene", seq_along(mock_stats))

  # Mock pathway (subset of gene names)
  mock_pathway <- names(mock_stats)[1:20]

  # Test with positive color value
  positive_color <- 1
  expect_no_error({
    positive_plot <- .plotEnrichment(
      pathway = mock_pathway,
      stats = mock_stats,
      color = positive_color,
      gseaParam = 1,
      ticksSize = 0.2
    )
  })

  # infinie stats
  expect_error({
    positive_plot <- .plotEnrichment(
      pathway = mock_pathway,
      stats = rep(-Inf, 200),
      color = positive_color,
      gseaParam = 1,
      ticksSize = 0.2
    )
  })

  # Verify plot is a ggplot object
  expect_true(inherits(positive_plot, "gg"))

  # Check layers
  expect_equal(length(positive_plot$layers), 1)

  # Verify layer is a geom_segment
  expect_true(
    inherits(positive_plot$layers[[1]]$geom, "GeomSegment")
  )

  # Test with negative color value
  negative_color <- -1
  expect_no_error({
    negative_plot <- .plotEnrichment(
      pathway = mock_pathway,
      stats = mock_stats,
      color = negative_color,
      gseaParam = 1,
      ticksSize = 0.2
    )
  })

  # Verify plot is a ggplot object
  expect_true(inherits(negative_plot, "gg"))
})
