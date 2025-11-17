# Test file for .barplot function
# Replace the content of tests/testthat/test-barplot.R with this content

library(testthat)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)

test_that(".barplot selects topN positives and negatives and assigns Type correctly", {
  RES_NES_strings <- data.frame(
    pathway = c("path1", "path2", "path3", "path4"),
    NES = c(3, -2, 1.5, -0.5),
    score = c(0.01, 0.05, 0.2, 0.9),
    stringsAsFactors = FALSE
  )

  p <- .barplot("dummy", type.sig = "score", topN = 1, RES_NES_strings = RES_NES_strings)
  expect_s3_class(p, "ggplot")

  # ggpubr stores the original data in p$data
  expect_true("Type" %in% names(p$data))
  expect_true("NES" %in% names(p$data))
  expect_true("pathway" %in% names(p$data))

  # Top1 positive is path1, top1 negative is path2
  expect_equal(sort(unique(p$data$Type)), c("protect", "risk"))
  expect_equal(p$data$Type[p$data$pathway == "path1"], "risk")
  expect_equal(p$data$Type[p$data$pathway == "path2"], "protect")

  # Should only have two rows (one positive, one negative)
  expect_equal(nrow(p$data), 2)

  # The wrapped label (y) should exist and correspond to pathway factor
  expect_true("y" %in% names(p$data))
  expect_true(all(p$data$pathway %in% c("path1", "path2")))
})

test_that(".barplot handles fewer than topN entries gracefully", {
  RES_NES_strings <- data.frame(
    pathway = c("p1", "p2"),
    NES = c(5, -5),
    score = c(0.1, 0.2),
    stringsAsFactors = FALSE
  )

  p <- .barplot("dummy", type.sig = "score", topN = 10, RES_NES_strings = RES_NES_strings)
  expect_s3_class(p, "ggplot")
  # both sides should appear since there are only one positive and one negative
  expect_equal(nrow(p$data), 2)
})

test_that(".barplot factor levels reflect reversed order of pathways", {
  RES_NES_strings <- data.frame(
    pathway = c("a", "b", "c"),
    NES = c(1, 2, -1),
    score = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )
  # topN=2: should pick NES positive top2 (a, b) and negative top1 (c) => combined a,b,c, reversed -> c,b,a
  p <- .barplot("dummy", type.sig = "score", topN = 2, RES_NES_strings = RES_NES_strings)
  expect_s3_class(p, "ggplot")
  df <- p$data

  # Check that pathway is a factor
  expect_true(is.factor(df$pathway))
  expected_levels <- c("c", "b", "a")  # rev(c("a","b","c"))
  expect_equal(levels(df$pathway), expected_levels)
})

test_that(".barplot errors when required columns missing", {
  # Missing the score column that type.sig refers to
  RES_NES_strings <- data.frame(
    pathway = c("x"),
    NES = c(1),
    stringsAsFactors = FALSE
  )
  expect_error(
    .barplot("dummy", type.sig = "score", topN = 1, RES_NES_strings = RES_NES_strings)
  )
})
