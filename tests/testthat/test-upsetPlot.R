library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

data("demo_GSE181574")
GSE181574.sigfun <- SummarizedExperiment::SummarizedExperiment(
    assays=list(abundance=as.matrix(expr.data)),
    rowData=S4Vectors::DataFrame(mapping, row.names=mapping$ensg_id),
    colData=SIG_MAT)
test_that("upsetPlot function works correctly", {

    result <- upsetPlot(GSE181574.sigfun)
    expect_type(result, "list")
    expect_named(result, c("upsetPlot", "Table_upsetPlot"))
    expect_s3_class(result$upsetPlot, "ggplot")
    expect_s3_class(result$Table_upsetPlot, "data.frame")

    result_bar <- upsetPlot(GSE181574.sigfun, type='bar')
    expect_s3_class(result_bar$upsetPlot, "ggplot")

    result_numeric <- upsetPlot(GSE181574.sigfun, showCategory=5)
    expect_s3_class(result_numeric$upsetPlot, "ggplot")

    specific_pathway <- c("REACTOME_CELL_CYCLE_MITOTIC",
                          "GOBP_CHROMOSOME_SEPARATION", "REACTOME_SYNTHESIS_OF_DNA")
    result_specific <- upsetPlot(GSE181574.sigfun, showCategory=specific_pathway)
    expect_s3_class(result_specific$upsetPlot, "ggplot")

})

test_that("upsetPlot handles errors correctly", {
    expect_snapshot(
        error = TRUE,
        upsetPlot(GSE181574.sigfun, type='invalid_method')
    )
})
