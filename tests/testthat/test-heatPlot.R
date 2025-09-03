library(testthat)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

data("demo_GSE181574")
GSE181574.sigfun <- SummarizedExperiment::SummarizedExperiment(
    assays=list(abundance=as.matrix(expr.data)),
    rowData=S4Vectors::DataFrame(mapping, row.names=mapping$ensg_id),
    colData=SIG_MAT)
test_that("heatPlot function works correctly", {

    result <- heatPlot(GSE181574.sigfun)
    expect_type(result, "list")
    expect_named(result, c("heatPlot", "table_heatPlot"))
    expect_s3_class(result$heatPlot, "ggplot")
    expect_s3_class(result$table_heatPlot, "data.frame")
    ### test color scale
    result_pos <- heatPlot(GSE181574.sigfun,showCategory=c('REACTOME_CELL_CYCLE_MITOTIC','REACTOME_DNA_REPLICATION'))
    result_all <- heatPlot(GSE181574.sigfun,showCategory=c('WP_RETINOBLASTOMA_GENE_IN_CANCER',
                                                           'REACTOME_CELL_CYCLE_MITOTIC','REACTOME_DNA_REPLICATION'))

    result_numeric <- heatPlot(GSE181574.sigfun, showCategory=5)
    expect_s3_class(result_numeric$heatPlot, "ggplot")

    specific_pathway <- c("REACTOME_CELL_CYCLE_MITOTIC",
                          "GOBP_CHROMOSOME_SEPARATION", "REACTOME_SYNTHESIS_OF_DNA")
    result_specific <- heatPlot(GSE181574.sigfun, showCategory=specific_pathway)
    expect_s3_class(result_specific$heatPlot, "ggplot")
})

test_that("heatPlot handles errors correctly", {
    expect_snapshot(
        error = TRUE,
        heatPlot(GSE181574.sigfun, dotSize='invalid_method')
    )
})
