# Test file for .upsetBarplot function
# Replace the content of tests/testthat/test-upsetBarplot.R with this content

library(testthat)
library(ggplot2)

# Ensure the SigFun namespace is loaded (without attaching)
if (!isNamespaceLoaded("SigFun")) requireNamespace("SigFun", quietly = TRUE)

# Create sample data for testing
sample_df <- data.frame(
  Path_Combination = c("A","A","B","B"),
  Count            = c(1,1,2,2),
  stringsAsFactors = FALSE
)

test_that(".upsetBarplot returns a ggplot object with geom_col", {
  p <- .upsetBarplot(sample_df)
  expect_s3_class(p, "ggplot")
  # First layer should be column geom
  expect_true(inherits(p$layers[[1]]$geom, "GeomCol"))
})

test_that(".upsetBarplot applies upsetTheme panel styling", {
  p <- .upsetBarplot(sample_df)
  bg <- p$theme$panel.background
  expect_s3_class(bg, "element_rect")
  expect_equal(bg$fill, "white")
  expect_equal(bg$colour, "black")
})

test_that(".upsetBarplot sets axis.text.y with size 12 and black color", {
  p <- .upsetBarplot(sample_df)
  ytext <- p$theme$axis.text.y
  expect_s3_class(ytext, "element_text")
  expect_equal(ytext$colour, "black")
  expect_equal(ytext$size, 12)
})

test_that(".upsetBarplot preserves blank x and y axis labels", {
  p <- .upsetBarplot(sample_df)
  expect_equal(p$labels$x, "")
  expect_equal(p$labels$y, "")
})
