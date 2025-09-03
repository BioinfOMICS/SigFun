library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

data("demo_GSE181574")
GSE181574.sigfun <- SummarizedExperiment::SummarizedExperiment(
    assays=list(abundance=as.matrix(expr.data)),
    rowData=S4Vectors::DataFrame(mapping, row.names=mapping$ensg_id),
    colData=SIG_MAT)
test_that("cnetPlot function works correctly", {

    result <- cnetPlot(GSE181574.sigfun)
    expect_type(result, "list")
    expect_named(result, c("cnetPlot", "TableEdge_cnetPlot", "TableNode_cnetPlot"))
    expect_s3_class(result$TableEdge_cnetPlot, "data.frame")

    result_color <- cnetPlot(GSE181574.sigfun, color=FALSE)
    expect_type(result_color, "list")
    expect_named(result_color, c("cnetPlot", "TableEdge_cnetPlot", "TableNode_cnetPlot"))
    expect_s3_class(result_color$TableEdge_cnetPlot, "data.frame")


    result_numeric <- cnetPlot(GSE181574.sigfun, showCategory=5)

    specific_pathway <- c("REACTOME_CELL_CYCLE_MITOTIC",
        "GOBP_CHROMOSOME_SEPARATION", "REACTOME_SYNTHESIS_OF_DNA")
    result_specific <- cnetPlot(GSE181574.sigfun, showCategory=specific_pathway)

})

test_that("cnetPlot handles errors correctly", {
  expect_snapshot(
    error = TRUE,
    cnetPlot(GSE181574.sigfun, dotColor='invalid_method')
  )
})

test_that("cnetPlot handles errors correctly", {
  expect_snapshot(
    error = TRUE,
    cnetPlot(GSE181574.sigfun, showCategory=c('invalid_pathway1', 'invalid_pathway2'))
  )
})
