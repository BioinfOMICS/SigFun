library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

data("demo_GSE181574")
GSE181574.sigfun <- SummarizedExperiment::SummarizedExperiment(
    assays=list(abundance=as.matrix(expr.data)),
    rowData=S4Vectors::DataFrame(mapping, row.names=mapping$ensg_id),
    colData=SIG_MAT)
test_that("gseaPlot function works correctly", {

    result <- gseaPlot(GSE181574.sigfun)
    expect_type(result, "list")


    specific_pathway <- c("REACTOME_CELL_CYCLE_MITOTIC",
          "GOBP_CHROMOSOME_SEPARATION", "REACTOME_SYNTHESIS_OF_DNA")
    result_specific <- gseaPlot(GSE181574.sigfun, showCategory=specific_pathway)

})
