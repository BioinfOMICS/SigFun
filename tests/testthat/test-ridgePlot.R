library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

data("demo_GSE181574")
GSE181574.sigfun <- SummarizedExperiment::SummarizedExperiment(
    assays=list(abundance=as.matrix(expr.data)),
    rowData=S4Vectors::DataFrame(mapping, row.names=mapping$ensg_id),
    colData=SIG_MAT)
test_that("ridgePlot function works correctly", {
    result <- ridgePlot(GSE181574.sigfun)
    expect_type(result, "list")
    expect_named(result, c("ridgePlot", "Table_ridgePlot"))
    expect_s3_class(result$ridgePlot, "ggplot")
    expect_s3_class(result$Table_ridgePlot, "data.frame")

    temp_file <- tempfile(fileext = ".png")
    suppressMessages(ggsave(temp_file, result$ridgePlot, width = 6, height = 4))
    file.remove(temp_file)

    result_numeric <- ridgePlot(GSE181574.sigfun, showCategory=5)
    expect_s3_class(result_numeric$ridgePlot, "ggplot")

    specific_pathway <- c("REACTOME_CELL_CYCLE_MITOTIC",
                          "GOBP_CHROMOSOME_SEPARATION", "REACTOME_SYNTHESIS_OF_DNA")
    result_specific <- ridgePlot(GSE181574.sigfun, showCategory=specific_pathway)
    expect_s3_class(result_specific$ridgePlot, "ggplot")

})
