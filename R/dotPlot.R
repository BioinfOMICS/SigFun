#' @title dotPlot
#' @description This function creates a dot plot visualization for Gene Set
#' Enrichment Analysis (GSEA) results stored in a SummarizedExperiment object.
#' The plot displays enriched gene sets with customizable visual encoding for
#' statistical significance and enrichment scores.
#' @param SE_data.fgsea A SummarizedExperiment object containing GSEA results.
#' The GSEA results should be stored in the metadata slot as 'gseaResult'.
#' @param showCategory An integer specifying the number of top enriched
#' categories to display in the plot, or a character vector of pathway names to
#' display specific pathways. Default is 10.
#' @param color A character vector specifying which variable to use for color
#' encoding. Must be one of "pvalue", "p.adjust", or "NES" (Normalized
#' Enrichment Score). Default is c('pvalue','p.adjust','NES').
#' @param breaklineN An integer specifying the number of characters after
#' which to break pathway names into multiple lines for better readability on
#' the y-axis. Default is \code{20}.
#' @return A list containing two elements:
#'   \describe{
#'     \item{\code{dotPlot}}{A ggplot2 object representing the dot plot with
#'           GSEA results. The plot includes the title "dotplot for GSEA" and
#'           shows the selected categories with dots sized by gene count and
#'           colored by the specified statistical measure (pvalue, p.adjust,
#'           or NES)}
#'     \item{\code{Table_dotPlot}}{A tibble/data.frame containing the plot data
#'           with columns: \code{original_name} (original pathway IDs),
#'           \code{label_name} (formatted pathway descriptions with line
#'           breaks),
#'           \code{Count} (number of genes in each pathway), \code{GeneRatio}
#'           (ratio of genes in pathway), and a column named after the
#'           \code{color} parameter containing the values used for color
#'           encoding}
#'   }
#' @export
#' @examples
#' data("demo_GSE181574")
#' dotPlot(SE_data.fgsea=GSE181574.sigfun)
dotPlot <- function(
        SE_data.fgsea, showCategory=10, color='pvalue', breaklineN=20){
    match.arg(color, c('pvalue', 'p.adjust', 'NES'))
    gseaRaw <- .extractDF(SE_data.fgsea, type='gseaRaw')
    gseaRaw@result <- gseaRaw@result |>
        dplyr::mutate(Description=.labelBreak(Description, breaklineN))
    dotPlot <- enrichplot::dotplot(
        gseaRaw, showCategory=showCategory, color=color) +
        ggplot2::ggtitle("dotplot for GSEA")
    Table_dotPlot <- dotPlot$data |>
        dplyr::select(original_name=ID, label_name=Description, Count,
            GeneRatio, dplyr::all_of(color)) |>
        tibble::remove_rownames()
    return(list(dotPlot=dotPlot, Table_dotPlot=Table_dotPlot))
}
