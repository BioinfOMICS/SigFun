# Test file for .upsetTheme function
# Replace the content of tests/testthat/test-upsetTheme.R with this content

library(testthat)
library(ggplot2)
library(grid)
# Ensure the SigFun namespace is loaded (without attaching)
if (!isNamespaceLoaded("SigFun")) requireNamespace("SigFun", quietly = TRUE)

test_that(".upsetTheme is a ggplot2 theme object", {
  th <- .upsetTheme
  expect_s3_class(th, "theme")
})

test_that("panel.background has white fill and black border", {
  bg <- .upsetTheme$panel.background
  expect_s3_class(bg, "element_rect")
  expect_equal(bg$fill, "white")
  expect_equal(bg$colour, "black")
})

test_that("axis.text.x and axis.ticks.y are blank elements", {
  th <- .upsetTheme
  expect_true(inherits(th$axis.text.x, "element_blank"))
  expect_true(inherits(th$axis.ticks.y, "element_blank"))
})

test_that("axis.text.y is element_text with black colour", {
  ytext <- .upsetTheme$axis.text.y
  expect_s3_class(ytext, "element_text")
  expect_equal(ytext$colour, "black")
})

test_that("axis.ticks.length is 0pt", {
  ttl <- .upsetTheme$axis.ticks.length
  expect_true(inherits(ttl, "unit"))
  expect_equal(as.numeric(ttl), 0)
  # as.character returns "0points" for zero-length unit
  expect_equal(as.character(ttl), "0points")
})

test_that("axis titles, axis.line, and panel.border are blank", {
  th <- .upsetTheme
  expect_true(inherits(th$axis.title.y, "element_blank"))
  expect_true(inherits(th$axis.title.x, "element_blank"))
  expect_true(inherits(th$axis.line, "element_blank"))
  expect_true(inherits(th$panel.border, "element_blank"))
})
