# Updated test file for .upsetTheme function (no border version)

library(testthat)
library(ggplot2)
library(grid)

# Ensure the SigFun namespace is loaded (without attaching)
if (!isNamespaceLoaded("SigFun")) requireNamespace("SigFun", quietly = TRUE)

test_that(".upsetTheme is a ggplot2 theme object", {
    th <- .upsetTheme
    expect_s3_class(th, "theme")
})

test_that("panel.background has white fill and no border color", {
    bg <- .upsetTheme$panel.background
    expect_s3_class(bg, "element_rect")
    expect_equal(bg$fill, "white")
    expect_true(is.null(bg$colour) || is.na(bg$colour))
})

test_that("axis.text.x is blank and axis.text.y is black text", {
    th <- .upsetTheme
    expect_true(inherits(th$axis.text.x, "element_blank"))
    ytext <- th$axis.text.y
    expect_s3_class(ytext, "element_text")
    expect_equal(ytext$colour, "black")
})

test_that("axis titles are blank elements", {
    th <- .upsetTheme
    expect_true(inherits(th$axis.title.x, "element_blank"))
    expect_true(inherits(th$axis.title.y, "element_blank"))
})

test_that("panel.border is blank", {
    th <- .upsetTheme
    expect_true(inherits(th$panel.border, "element_blank"))
})

# Optional: relaxed check on axis.line, only test if exists
test_that("axis.line is either blank or not set", {
    th <- .upsetTheme
    if (!is.null(th$axis.line)) {
        expect_true(inherits(th$axis.line, "element_blank"))
    } else {
        succeed()
    }
})
