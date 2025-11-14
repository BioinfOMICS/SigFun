# Test file for .upsetBoxplot function
# Updated for new conditional fillColor handling

library(testthat)
library(ggplot2)

# Ensure the SigFun namespace is loaded (without attaching)
if (!isNamespaceLoaded("SigFun")) requireNamespace("SigFun", quietly = TRUE)

# Create sample data for testing
set.seed(123)
sample_df <- data.frame(
    Path_Combination = rep(c("A", "B"), each = 5),
    Coef             = c(rnorm(5), rnorm(5, mean = 2)),
    stringsAsFactors = FALSE
)

test_that(".upsetBoxplot returns a ggplot object with boxplot and jitter layers", {
    p <- .upsetBoxplot(sample_df)
    expect_s3_class(p, "ggplot")
    # First layer: boxplot
    expect_true(inherits(p$layers[[1]]$geom, "GeomBoxplot"))
    # Second layer: jittered points
    expect_true(inherits(p$layers[[2]]$geom, "GeomPoint"))
    expect_true(inherits(p$layers[[2]]$position, "PositionJitter"))
})

test_that(".upsetBoxplot works without fillColor argument (no warning)", {
    expect_silent(.upsetBoxplot(sample_df))
})

test_that(".upsetBoxplot applies custom fillColor correctly", {
    p <- .upsetBoxplot(sample_df, fillColor = "#5979A3")
    expect_equal(p$layers[[1]]$aes_params$fill, "#5979A3")
})

test_that(".upsetBoxplot applies upsetTheme background styling", {
    p <- .upsetBoxplot(sample_df)
    bg <- p$theme$panel.background
    expect_s3_class(bg, "element_rect")
    expect_equal(bg$fill, "white")
    expect_true(is.null(bg$colour) || is.na(bg$colour) || bg$colour == "black")
})

test_that(".upsetBoxplot sets axis.text.y style correctly", {
    p <- .upsetBoxplot(sample_df, fontSize = 12)
    ytext <- p$theme$axis.text.y
    expect_s3_class(ytext, "element_text")
    expect_equal(ytext$colour, "black")
    expect_equal(ytext$size, 12)
})

test_that(".upsetBoxplot keeps x and y axis labels blank", {
    p <- .upsetBoxplot(sample_df)
    expect_equal(p$labels$x, "")
    expect_equal(p$labels$y, "")
})
