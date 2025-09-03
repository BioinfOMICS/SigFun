library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

data("demo_GSE181574")
GSE181574.sigfun <- SummarizedExperiment::SummarizedExperiment(
    assays=list(abundance=as.matrix(expr.data)),
    rowData=S4Vectors::DataFrame(mapping, row.names=mapping$ensg_id),
    colData=SIG_MAT)

test_that("emapPlot function works correctly", {

    result <- emapPlot(GSE181574.sigfun)
    expect_type(result, "list")
    expect_named(result, c("emapPlot", "Table_emapPlot"))
    expect_s3_class(result$emapPlot, "ggplot")
    expect_s3_class(result$Table_emapPlot, "data.frame")

    result_numeric <- emapPlot(GSE181574.sigfun, showCategory=5)
    expect_named(result_numeric, c("emapPlot", "Table_emapPlot"))
    expect_s3_class(result_numeric$emapPlot, "ggplot")
    expect_s3_class(result_numeric$Table_emapPlot, "data.frame")

    specific_pathway <- c("REACTOME_CELL_CYCLE_MITOTIC",
          "GOBP_CHROMOSOME_SEPARATION", "REACTOME_SYNTHESIS_OF_DNA")
    result_specific <- emapPlot(GSE181574.sigfun, showCategory=specific_pathway)
    expect_named(result_specific, c("emapPlot", "Table_emapPlot"))
    expect_s3_class(result_specific$emapPlot, "ggplot")
    expect_s3_class(result_specific$Table_emapPlot, "data.frame")

})
