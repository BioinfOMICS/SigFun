# Test file for .cowplotText function
# Replace the content of tests/testthat/test-cowplotText.R with this content

library(testthat)
library(cowplot)
library(grid)  # for is.grob

# Ensure the package is properly loaded (if .cowplotText is in SigFun namespace)
if (!isNamespaceLoaded("SigFun")) {
  requireNamespace("SigFun", quietly = TRUE)
}

# helper to recursively extract text labels from a grob
extract_text_labels <- function(grob) {
  labels <- character(0)
  if (inherits(grob, "text")) {
    labels <- c(labels, as.character(grob$label))
  }
  if (!is.null(grob$children)) {
    for (child in grob$children) {
      labels <- c(labels, extract_text_labels(child))
    }
  }
  labels
}

get_all_text_labels <- function(plot_obj) {
  g <- ggplot2::ggplotGrob(plot_obj)
  unique(unlist(lapply(g$grobs, extract_text_labels)))
}

test_that(".cowplotText basic functionality returns a ggplot object", {
  p <- .cowplotText("Hello", list())
  expect_s3_class(p, "gg")        # gg-compatible
  expect_s3_class(p, "ggplot")
  # Should contain the string "Hello" somewhere in its text grobs
  labels <- get_all_text_labels(p)
  expect_true(any(grepl("Hello", labels, fixed = TRUE)))
})

test_that(".cowplotText accepts custom style parameters and includes text grob", {
  style <- list(size = 30, colour = "red", fontface = "bold", x = 0.3, y = 0.7, hjust = 0, vjust = 1)
  p <- .cowplotText("Custom Style", style)
  expect_s3_class(p, "ggplot")
  labels <- get_all_text_labels(p)
  expect_true(any(grepl("Custom Style", labels, fixed = TRUE)))
})

test_that(".cowplotText handles long ASCII strings without error", {
  long_text <- "This is a very long ASCII string with special characters !@#$%^&*()-_=+[]{};:'\",.<>/? and numbers 1234567890 to check robustness."
  p <- .cowplotText(long_text, list(size = 15))
  expect_s3_class(p, "ggplot")
  labels <- get_all_text_labels(p)
  # match a substring that actually exists
  expect_true(any(grepl("This is a very long ASCII string", labels, fixed = TRUE)))
})

test_that(".cowplotText coerces numeric input to string", {
  p_num <- .cowplotText(12345, list())
  expect_s3_class(p_num, "ggplot")
  labels <- get_all_text_labels(p_num)
  expect_true(any(grepl("12345", labels, fixed = TRUE)))
})

test_that(".cowplotText output differs when input text differs", {
  p1 <- .cowplotText("A", list(size = 20))
  p2 <- .cowplotText("B", list(size = 20))
  labels1 <- get_all_text_labels(p1)
  labels2 <- get_all_text_labels(p2)
  expect_false(identical(labels1, labels2))
})

test_that(".cowplotText invalid non-numeric size produces error", {
  expect_error(.cowplotText("Bad size", list(size = "big")))
})

test_that(".cowplotText handles NA size gracefully", {
  p <- .cowplotText("NA size", list(size = NA))
  expect_s3_class(p, "ggplot")
  labels <- get_all_text_labels(p)
  expect_true(any(grepl("NA size", labels, fixed = TRUE)))
})

