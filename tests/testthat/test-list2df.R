# Test file for .list2df function
# Replace the content of tests/testthat/test-list2df.R with this content

library(testthat)
# Ensure SigFun namespace is loaded
if (!isNamespaceLoaded("SigFun"))
  requireNamespace("SigFun", quietly = TRUE)

test_that(".list2df converts a simple named list to a data.frame correctly", {
  inputList <- list(
    A = c("gA1", "gA2"),
    B = c("gB")
  )
  out <- .list2df(inputList)
  expected <- data.frame(
    categoryID    = c("A", "A", "B"),
    Gene          = c("gA1", "gA2", "gB"),
    stringsAsFactors = FALSE
  )
  rownames(expected) <- NULL
  expect_equal(out, expected)
})

test_that(".list2df handles an empty input list by returning NULL", {
  inputList <- list()
  out <- .list2df(inputList)
  expect_null(out)
})

test_that(".list2df skips list elements of length zero", {
  inputList <- list(
    A = character(0),
    B = c("g1")
  )
  out <- .list2df(inputList)
  expected <- data.frame(
    categoryID    = c("B"),
    Gene          = c("g1"),
    stringsAsFactors = FALSE
  )
  rownames(expected) <- NULL
  expect_equal(out, expected)
})

test_that(".list2df works with numeric gene vectors", {
  inputList <- list(
    X = 1:3
  )
  out <- .list2df(inputList)
  expected <- data.frame(
    categoryID    = rep("X", 3),
    Gene          = 1:3,
    stringsAsFactors = FALSE
  )
  rownames(expected) <- NULL
  expect_equal(out, expected)
})
