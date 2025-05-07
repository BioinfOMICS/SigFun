library(testthat)
library(ggplot2)

test_that(".nesbarplot function works correctly", {
  # Test positive NES value
  positive_nes <- 0.75
  positive_maxNES <- 1

  positive_plot <- .nesbarplot(positive_nes, positive_maxNES)

  # Check that the plot is a ggplot object
  expect_true(inherits(positive_plot, "gg"))

  # Verify color for positive NES (should be reddish)
  expect_true(all(
    positive_plot$layers[[1]]$aes_params$fill == '#FD7C69'
  ))

  # Verify x-axis limits
  x_limits <- layer_scales(positive_plot)$x$limits
  expect_equal(x_limits, c(0, positive_maxNES))

  # Test negative NES value
  negative_nes <- -0.5
  negative_maxNES <- 1

  negative_plot <- .nesbarplot(negative_nes, negative_maxNES)

  # Check that the plot is a ggplot object
  expect_true(inherits(negative_plot, "gg"))

  # Verify color for negative NES (should be bluish)
  expect_true(all(
    negative_plot$layers[[1]]$aes_params$fill == '#5F90BB'
  ))

  # Verify x-axis limits for negative case
  x_limits_neg <- layer_scales(negative_plot)$x$limits
  expect_equal(x_limits_neg, c(0, negative_maxNES))

  # Check theme elements are blank
  theme_elements <- c(
    "panel.background",
    "plot.background",
    "axis.text.x",
    "axis.line",
    "axis.text",
    "axis.ticks",
    "panel.grid",
    "axis.title"
  )

  for (element in theme_elements) {
    expect_true(
      inherits(
        positive_plot$theme[[element]],
        "element_blank"
      )
    )
  }
})
