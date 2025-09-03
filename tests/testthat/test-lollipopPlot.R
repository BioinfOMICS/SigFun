library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

data("demo_GSE181574")
GSE181574.sigfun <- SummarizedExperiment::SummarizedExperiment(
    assays=list(abundance=as.matrix(expr.data)),
    rowData=S4Vectors::DataFrame(mapping, row.names=mapping$ensg_id),
    colData=SIG_MAT)
test_that("lollipopPlot function works correctly", {

    result <- lollipopPlot(GSE181574.sigfun)
    expect_type(result, "list")
    expect_named(result, c("lollipopPlot", "Table_lollipopPlot"))
    expect_s3_class(result$lollipopPlot, "ggplot")
    expect_s3_class(result$Table_lollipopPlot, "data.frame")

    temp_file <- tempfile(fileext = ".png")
    suppressMessages(ggsave(temp_file, result$lollipopPlot, width = 6, height = 4))
    file.remove(temp_file)


    result_numeric <- lollipopPlot(GSE181574.sigfun, showCategory=5)
    expect_s3_class(result_numeric$lollipopPlot, "ggplot")
    expect_lte(nrow(result_numeric$Table_lollipopPlot), 5)

    specific_pathway <- c("REACTOME_CELL_CYCLE_MITOTIC",
                          "GOBP_CHROMOSOME_SEPARATION", "REACTOME_SYNTHESIS_OF_DNA")
    result_specific <- lollipopPlot(GSE181574.sigfun, showCategory=specific_pathway)
    expect_s3_class(result_specific$lollipopPlot, "ggplot")
    expect_lte(nrow(result_specific$Table_lollipopPlot),
        length(specific_pathway))

})
