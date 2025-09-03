# Test file for .scale_bar function
# Replace the content of tests/testthat/test-scale_bar.R with this content

library(testthat)
library(ggplot2)

# Ensure the SigFun namespace is loaded (without attaching)
if (!isNamespaceLoaded("SigFun")) requireNamespace("SigFun", quietly = TRUE)

test_that(".scale_bar returns scale_fill_gradient for non-negative values", {
  s1 <- .scale_bar(c(0, 1, 5))
  call1 <- s1$call
  # Extract function name (last part of the call)
  fn1 <- tail(as.character(call1[[1]]), 1)
  expect_equal(fn1, "scale_fill_gradient")
  expect_equal(call1$low, "#FF5151")
  expect_equal(call1$high, "#F8806A")
})

test_that(".scale_bar returns scale_fill_gradient for non-positive values", {
  s2 <- .scale_bar(c(-5, -1, 0))
  call2 <- s2$call
  fn2 <- tail(as.character(call2[[1]]), 1)
  expect_equal(fn2, "scale_fill_gradient")
  expect_equal(call2$low, "#4169E1")
  expect_equal(call2$high, "#deebf7")
})

test_that(".scale_bar returns scale_fill_gradient2 for values spanning zero", {
  s3 <- .scale_bar(c(-3, 0, 4))
  call3 <- s3$call
  fn3 <- tail(as.character(call3[[1]]), 1)
  expect_equal(fn3, "scale_fill_gradient2")
  expect_equal(call3$low, "#4169E1")
  expect_equal(call3$mid, "white")
  expect_equal(call3$high, "#FF5151")
})
