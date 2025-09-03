# Test file for .upsetBoxplot function
# Replace the content of tests/testthat/test-upsetBoxplot.R with this content

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

test_that(".upsetBoxplot returns a ggplot object with boxplot and jitter", {
  p <- .upsetBoxplot(sample_df)
  expect_s3_class(p, "ggplot")

  # First layer is boxplot
  expect_true(inherits(p$layers[[1]]$geom, "GeomBoxplot"))
  # Second layer is point geom with jitter position
  expect_true(inherits(p$layers[[2]]$geom, "GeomPoint"))
  expect_true(inherits(p$layers[[2]]$position, "PositionJitter"))
})

test_that(".upsetBoxplot retains upsetTheme panel background styling", {
  p <- .upsetBoxplot(sample_df)
  bg <- p$theme$panel.background
  expect_s3_class(bg, "element_rect")
  expect_equal(bg$fill, "white")
  expect_equal(bg$colour, "black")
})

test_that(".upsetBoxplot sets axis.text.y with size 12 and black color", {
  p <- .upsetBoxplot(sample_df)
  ytext <- p$theme$axis.text.y
  expect_s3_class(ytext, "element_text")
  expect_equal(ytext$colour, "black")
  expect_equal(ytext$size, 12)
})

test_that(".upsetBoxplot preserves blank x and y axis labels", {
  p <- .upsetBoxplot(sample_df)
  expect_equal(p$labels$x, "")
  expect_equal(p$labels$y, "")
})
