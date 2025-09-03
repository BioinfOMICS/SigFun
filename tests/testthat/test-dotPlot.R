library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

data("demo_GSE181574")
GSE181574.sigfun <- SummarizedExperiment::SummarizedExperiment(
    assays=list(abundance=as.matrix(expr.data)),
    rowData=S4Vectors::DataFrame(mapping, row.names=mapping$ensg_id),
    colData=SIG_MAT)

test_that("dotPlot function works correctly", {

    result <- dotPlot(GSE181574.sigfun)
    expect_type(result, "list")
    expect_named(result, c("dotPlot", "Table_dotPlot"))
    expect_s3_class(result$Table_dotPlot, "data.frame")

    result_numeric <- dotPlot(GSE181574.sigfun, showCategory=5)
    expect_named(result_numeric, c("dotPlot", "Table_dotPlot"))
    expect_s3_class(result_numeric$Table_dotPlot, "data.frame")

    specific_pathway <- c("REACTOME_CELL_CYCLE_MITOTIC",
          "GOBP_CHROMOSOME_SEPARATION", "REACTOME_SYNTHESIS_OF_DNA")
    result_specific <- dotPlot(GSE181574.sigfun, showCategory=specific_pathway)
    expect_named(result_specific, c("dotPlot", "Table_dotPlot"))
    expect_s3_class(result_specific$Table_dotPlot, "data.frame")

})

test_that("dotPlot handles errors correctly", {
    expect_snapshot(
        error = TRUE,
        dotPlot(GSE181574.sigfun, color='invalid_method')
    )
})
