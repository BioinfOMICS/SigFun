library(testthat)
library(ggplot2)

# Ensure the package is properly loaded
if (!isNamespaceLoaded("SigFun")) {
    requireNamespace("SigFun", quietly = TRUE)
}

test_that(".nesbarplot function works correctly", {
    # ---- Positive NES ----
    positive_nes <- 0.75
    positive_maxNES <- 1
    positive_plot <- .nesbarplot(positive_nes, positive_maxNES)

    # Check ggplot object
    expect_s3_class(positive_plot, "ggplot")

    # Verify fill color for positive NES
    expect_true(all(positive_plot$layers[[1]]$aes_params$fill == '#D25C43'))

    # Verify x-axis limits
    x_limits <- layer_scales(positive_plot)$x$limits
    expect_equal(x_limits, c(0, positive_maxNES))

    # ---- Negative NES ----
    negative_nes <- -0.5
    negative_maxNES <- 1
    negative_plot <- .nesbarplot(negative_nes, negative_maxNES)

    # Check ggplot object
    expect_s3_class(negative_plot, "ggplot")

    # Verify fill color for negative NES
    expect_true(all(negative_plot$layers[[1]]$aes_params$fill == '#5979A3'))

    # Verify x-axis limits
    x_limits_neg <- layer_scales(negative_plot)$x$limits
    expect_equal(x_limits_neg, c(0, negative_maxNES))

    # ---- Theme Checks ----
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
        expect_true(inherits(positive_plot$theme[[element]], "element_blank"))
    }
})
