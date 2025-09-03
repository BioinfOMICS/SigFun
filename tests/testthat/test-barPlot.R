library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

data("demo_GSE181574")
GSE181574.sigfun <- SummarizedExperiment::SummarizedExperiment(
    assays=list(abundance=as.matrix(expr.data)),
    rowData=S4Vectors::DataFrame(mapping, row.names=mapping$ensg_id),
    colData=SIG_MAT)

test_that("barPlot function works correctly", {

    result <- barPlot(GSE181574.sigfun)
    expect_type(result, "list")
    expect_named(result, c("barPlot", "Table_barPlot"))
    expect_s3_class(result$barPlot, "ggplot")
    expect_s3_class(result$Table_barPlot, "data.frame")

    for (color_param in c("pvalue", "p.adjust", "qvalue")) {
        result_color <- barPlot(GSE181574.sigfun, alpha=color_param)
        expect_s3_class(result_color$barPlot, "ggplot")
        expect_true(color_param %in% colnames(result_color$Table_barPlot))
    }
})

test_that("barPlot handles errors correctly", {
    expect_snapshot(
        error = TRUE,
        barPlot(GSE181574.sigfun, alpha='invalid_method')
    )
})

