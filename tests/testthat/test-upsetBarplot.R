# tests/testthat/test-upsetBarplot.R

library(testthat)
library(ggplot2)

if (!isNamespaceLoaded("SigFun")) requireNamespace("SigFun", quietly = TRUE)

sample_df <- data.frame(
    Path_Combination = c("A", "A", "B", "B"),
    Count            = c(1, 1, 2, 2),
    stringsAsFactors = FALSE
)

test_that(".upsetBarplot returns a ggplot object with geom_col", {
    p <- .upsetBarplot(sample_df, fontSize = 10, fillColor = "#5979A3")
    expect_s3_class(p, "ggplot")
    expect_true(inherits(p$layers[[1]]$geom, "GeomCol"))
})

test_that(".upsetBarplot uses provided fill color", {
    fill_color <- "#5979A3"
    p <- .upsetBarplot(sample_df, fontSize = 10, fillColor = fill_color)
    geom_params <- p$layers[[1]]$aes_params
    expect_equal(geom_params$fill, fill_color)
})

test_that(".upsetBarplot applies .upsetTheme and correct axis styling", {
    p <- .upsetBarplot(sample_df, fontSize = 12, fillColor = "#5979A3")
    ytext <- p$theme$axis.text.y
    expect_s3_class(ytext, "element_text")
    expect_equal(ytext$colour, "black")
    expect_equal(ytext$size, 12)
    expect_s3_class(p$theme$axis.ticks.x, "element_blank")
    axis_line <- p$theme$axis.line
    expect_s3_class(axis_line, "element_line")
    expect_equal(axis_line$colour, "grey40")
})

test_that(".upsetBarplot preserves empty x and y labels", {
    p <- .upsetBarplot(sample_df)
    expect_equal(p$labels$x, "")
    expect_equal(p$labels$y, "")
})

test_that(".upsetBarplot scales y-axis correctly", {
    p <- .upsetBarplot(sample_df)
    expect_true("ScaleContinuousPosition" %in% class(p$scales$scales[[1]]))
})
