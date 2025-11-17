# Test file for .upsetPlot function
# Replace the content of tests/testthat/test-upsetplot.R with this content
library(testthat)
library(ggplot2)

# Ensure the SigFun namespace is loaded (without attaching)
if (!isNamespaceLoaded("SigFun")) requireNamespace("SigFun", quietly = TRUE)

test_that(".upsetPlot returns a ggplot object with correct structure", {
  df <- data.frame(
    Path_Combination = c("A", "B"),
    P1 = c(TRUE, FALSE),
    P2 = c(FALSE, TRUE),
    P3 = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  p <- .upsetPlot(df)

  # Should be a ggplot
  expect_s3_class(p, "ggplot")

  # Check if the plot has the expected layers
  expect_true(length(p$layers) >= 2)  # Should have geom_point and geom_line

  # Check that a discrete color scale is present
  scale_color <- p$scales$get_scales("colour")
  if (is.null(scale_color)) scale_color <- p$scales$get_scales("color")
  expect_true(!is.null(scale_color))
  expect_true(inherits(scale_color, "ScaleDiscrete"))
})

test_that(".upsetPlot has correct color mapping on build", {
  df <- data.frame(
    Path_Combination = c("A", "B"),
    P1 = c(TRUE, FALSE),
    P2 = c(FALSE, TRUE),
    P3 = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  p <- .upsetPlot(df)
  built <- ggplot_build(p)

  # First layer data should have a color aesthetic
  d <- built$data[[1]]
  expect_true("colour" %in% names(d) || "color" %in% names(d))
})

test_that("axis.text.y size is 12 when fewer than 10 pathways", {
  df <- data.frame(
    Path_Combination = c("A", "B"),
    X = c(TRUE, FALSE),
    Y = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  p <- .upsetPlot(df)
  th_y <- p$theme$axis.text.y
  expect_s3_class(th_y, "element_text")
  expect_equal(th_y$size, 12)
})

test_that("axis.text.y size is 8 when 10 or more pathways", {
  mat <- as.data.frame(matrix(TRUE, nrow = 1, ncol = 10))
  names(mat) <- paste0("L", seq_len(10))
  df <- cbind(Path_Combination = "X", mat, stringsAsFactors = FALSE)

  p <- .upsetPlot(df)
  th_y <- p$theme$axis.text.y
  expect_s3_class(th_y, "element_text")
  expect_equal(th_y$size, 8)
})

test_that("Pathway factor levels are reversed unique order", {
  df <- data.frame(
    Path_Combination = c("A", "B"),
    Alpha = c(TRUE, FALSE),
    Beta  = c(FALSE, TRUE),
    Gamma = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  p <- .upsetPlot(df)
  lvl <- levels(p$data$Pathway)
  expect_equal(lvl, c("Gamma", "Beta", "Alpha"))
})

test_that(".upsetPlot transforms data correctly", {
  df <- data.frame(
    Path_Combination = c("A", "B"),
    P1 = c(TRUE, FALSE),
    P2 = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  p <- .upsetPlot(df)
  plot_data <- p$data
  expect_true(all(c("Path_Combination","Pathway","value") %in% names(plot_data)))
  expect_true(is.factor(plot_data$Pathway))
  expect_true(is.logical(plot_data$value))
})

test_that(".upsetPlot handles single pathway correctly", {
  df <- data.frame(
    Path_Combination = c("A", "B"),
    OnlyPath = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  p <- .upsetPlot(df)
  expect_s3_class(p, "ggplot")
  expect_equal(nlevels(p$data$Pathway), 1)
  expect_equal(levels(p$data$Pathway), "OnlyPath")
})
