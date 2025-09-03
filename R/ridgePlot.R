#' @title ridgePlot
#' @description The \code{ridgePlot} function creates a ridge plot visualization
#' for Gene Set Enrichment Analysis (GSEA) results. Ridge plots display the
#' distribution of gene expression values for each gene set as overlapping
#' density curves, providing an intuitive way to visualize and compare the
#' enrichment patterns across multiple pathways. Each ridge represents a gene
#' set, with the curve shape indicating the distribution of genes within that
#' pathway along the ranked gene list.
#' @param SE_data.fgsea A \code{SummarizedExperiment} object containing GSEA
#'   results in its metadata. The object must contain:
#'   \itemize{
#'     \item \code{metadata(SE_data.fgsea)$gseaResult}: GSEA result object
#'           compatible with enrichplot functions, containing pathway
#'           information, gene sets, and enrichment statistics
#'   }
#' @param showCategory An integer specifying the number of top-ranked gene sets
#' to display in the ridge plot. Gene sets are typically ordered by statistical
#' significance (p-value or adjusted p-value) or enrichment score.
#' Default is \code{10}.
#' @param breaklineN An integer specifying the number of characters after which
#' to break pathway names into multiple lines for better readability on the
#' y-axis. Default is \code{30}.
#' @return A list containing two elements:
#'   \describe{
#'     \item{\code{ridgePlot}}{A ggplot2 object containing the ridge plot
#'           visualization with formatted pathway labels}
#'     \item{\code{Table_ridgePlot}}{A tibble/data.frame containing the plot
#'           data with columns: \code{original_name} (original pathway names),
#'           \code{label_name} (formatted pathway names with line breaks),
#'           \code{p.adjust} (adjusted p-values), and \code{value} (enrichment
#'           values used for ridge distribution)}
#'   }
#' @export
#' @examples
#' data("demo_GSE181574")
#' ridgePlot(SE_data.fgsea=GSE181574.sigfun)
ridgePlot <- function(SE_data.fgsea, showCategory=10, breaklineN=30){
    gseaRaw <- .extractDF(SE_data.fgsea, type='gseaRaw')
    ridgePlot <- enrichplot::ridgeplot(gseaRaw, showCategory=showCategory)  +
        ggplot2::scale_y_discrete(labels=function(x){
            .labelBreak(x, breaklineN)})
    Table_ridgePlot <- ridgePlot$data |>
        dplyr::mutate(label_name=.labelBreak(category, breaklineN)) |>
        dplyr::select(original_name=category, label_name, p.adjust, value) |>
        tibble::remove_rownames()
    return(list(ridgePlot=ridgePlot, Table_ridgePlot=Table_ridgePlot))
}
