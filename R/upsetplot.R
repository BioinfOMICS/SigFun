#' @title upsetPlot
#' @description The \code{upsetPlot} function generates an UpSet plot for
#' visualizing overlapping genes across multiple pathways from gene set
#' enrichment analysis(GSEA) results. The plot consists of two panels: a boxplot
#' showing the distribution of correlation coefficients for each pathway
#' combination, and an intersection matrix showing which pathways contribute to
#' each combination.
#' @param SE_data.fgsea A \code{SummarizedExperiment} object containing GSEA
#'   results and correlation data. The object must contain:
#'   \itemize{
#'     \item \code{metadata(SE_data.fgsea)$gseaResult}: GSEA result object with
#'           pathway information and gene sets
#'   }
#' @param showCategory Either a numeric value indicating the number of top
#' pathways to display, or a character vector specifying particular pathway
#' names to show. Default is \code{5}.
#' @param breaklineN An integer specifying the number of characters after which
#' to break pathway names into multiple lines for better readability.
#' Default is \code{30}.
#' @param type A character string specifying the type of plot to display in the
#' upper panel. Must be one of \code{c('bar', 'box')}. Default is \code{'box'}.
#' @return A list containing two elements:
#'   \describe{
#'     \item{\code{upsetPlot}}{A ggplot2 object containing the UpSet plot
#'           visualization with two panels: the upper panel shows either
#'           bar plot or box plot (depending on \code{type} parameter)
#'           displaying pathway combination statistics, and the lower panel
#'           shows the intersection matrix indicating which pathways
#'           contribute to each combination}
#'     \item{\code{Table_upsetPlot}}{A data.frame containing the raw data
#'           used for generating the UpSet plot, including pathway combinations,
#'           gene overlaps, and associated statistics}
#'   }
#' @export
#' @examples
#' data("demo_GSE181574")
#' upsetPlot(SE_data.fgsea=GSE181574.sigfun)
upsetPlot <- function(SE_data.fgsea, showCategory=5, breaklineN=30, type='box'){
    match.arg(type, c('bar', 'box'))
    gseaRaw <- .extractDF(SE_data.fgsea, type="gseaRaw")
    geneSets <- .extractGeneSets(gseaRaw, showCategory)
    fontSize <- ifelse(length(geneSets) < 10, 12, 8)
    upsetData <- .upsetData(geneSets, gseaRaw, breaklineN, type)
    upsetTop <- switch(type,
        bar=.upsetBarplot(upsetData$data),
        box=.upsetBoxplot(upsetData$data))
    upsetdottom <-.upsetPlot(upsetData$data)
    n.path <- length(unique(geneSets$categoryID))
    if(n.path < 5){
        heights <- c(1-n.path/20, n.path/20)
    }else{
        heights <- c(0.3, 0.7)
    }
    upsetPlot <- ggpubr::ggarrange(
        upsetTop, upsetdottom, ncol=1, align="v", heights=heights)
    Table_upsetPlot <- upsetData$Rawdata
    return(list(upsetPlot=upsetPlot, Table_upsetPlot=Table_upsetPlot))
}
