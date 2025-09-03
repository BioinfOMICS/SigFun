# Test file for .scale_fill function
# Replace the content of tests/testthat/test-scale_fill.R with this content

library(testthat)
library(ggplot2)
library(scales)

# Ensure the SigFun namespace is loaded (without attaching)
if (!isNamespaceLoaded("SigFun")) requireNamespace("SigFun", quietly = TRUE)

test_that(".scale_fill for non-negative values uses zero->pos mapping", {
  vals <- c(0, 2, 5)
  s1   <- .scale_fill(vals)
  call1 <- s1$call

  # Function should be scale_fill_gradient2
  fn1 <- tail(as.character(call1[[1]]), 1)
  expect_equal(fn1, "scale_fill_gradient2")

  # low and high arguments refer to zero and pos parameters
  expect_equal(as.character(call1$low), "zero")
  expect_equal(as.character(call1$high), "pos")

  # limits evaluates to c(0, max(vals))
  lims1 <- eval(call1$limits, envir = list(value = vals))
  expect_equal(lims1, c(0, max(vals)))

  # oob should be scales::squish
  oob1_parts <- as.character(call1$oob)
  expect_equal(tail(oob1_parts, 1), "squish")
})

test_that(".scale_fill for non-positive values uses neg->zero mapping", {
  vals <- c(-5, -1, 0)
  s2   <- .scale_fill(vals)
  call2 <- s2$call

  fn2 <- tail(as.character(call2[[1]]), 1)
  expect_equal(fn2, "scale_fill_gradient2")

  expect_equal(as.character(call2$low), "neg")
  expect_equal(as.character(call2$high), "zero")

  lims2 <- eval(call2$limits, envir = list(value = vals))
  expect_equal(lims2, c(min(vals), 0))

  oob2_parts <- as.character(call2$oob)
  expect_equal(tail(oob2_parts, 1), "squish")
})

test_that(".scale_fill for values spanning zero uses neg->mid->pos mapping", {
  vals <- c(-3, 0, 4)
  s3   <- .scale_fill(vals)
  call3 <- s3$call

  fn3 <- tail(as.character(call3[[1]]), 1)
  expect_equal(fn3, "scale_fill_gradient2")

  expect_equal(as.character(call3$low), "neg")
  expect_equal(as.character(call3$mid), "zero")
  expect_equal(as.character(call3$high), "pos")

  lims3 <- eval(call3$limits, envir = list(value = vals))
  expect_equal(lims3, c(min(vals), max(vals)))

  oob3_parts <- as.character(call3$oob)
  expect_equal(tail(oob3_parts, 1), "squish")
})
