library(testthat)
library(cowplot)
library(ggplot2)

test_that(".cowplotText function works correctly", {
  # Test with basic text and style
  basic_style <- list(
    x=0.5,
    y=0.5,
    size=12,
    color="black"
  )

  plot <- .cowplotText("Test Text", basic_style)

  # Check that the plot is a ggplot object
  expect_true(inherits(plot, "gg"))

  # Check that the plot is created by ggdraw
  expect_true(inherits(plot$layers[[1]]$geom, "GeomText"))

  # Test with additional style parameters
  complex_style <- list(
    x=0.3,
    y=0.7,
    size=16,
    color="red",
    fontface="bold"
  )

  complex_plot <- .cowplotText("Complex Text", complex_style)

  # Verify complex plot attributes
  expect_true(inherits(complex_plot, "gg"))

  # Optional: Could add more specific checks about text properties if needed
})
