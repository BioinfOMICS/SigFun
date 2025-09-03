#' @title lollipopPlot
#' @description The \code{lollipopPlot} function creates a lollipop plot
#' (enhanced dot plot) to visualize gene set enrichment analysis (GSEA) results.
#'  This visualization displays pathways as horizontal bars with dots at the
#' end, where the position of each dot represents the Normalized Enrichment
#' Score (NES) and additional visual elements indicate statistical
#' significance. The plot provides an intuitive way to compare pathway
#' enrichment across multiple biological processes, with pathways ordered by
#' their enrichment scores.
#' @param SE_data.fgsea A \code{SummarizedExperiment} object containing GSEA
#'   results in its metadata. The object must contain:
#'   \itemize{
#'     \item \code{metadata(SE_data.fgsea)$gseaResult}: GSEA result object with
#'     pathway information, gene sets, enrichment statistics, and NES values
#'   }
#' @param showCategory Either an integer specifying the number of top pathways
#'   to display (ranked by significance), or a character vector of specific
#'   pathway IDs to include in the visualization. Default is \code{10}.
#' @param breaklineN An integer specifying the maximum number of characters per
#'   line before wrapping pathway names to the next line for better readability.
#'   Underscores in pathway names are automatically replaced with spaces.
#'   Default is \code{30}.
#' @return A list containing two elements:
#'   \describe{
#'     \item{\code{lollipopPlot}}{A ggplot2 object containing the lollipop plot
#'           visualization showing pathway enrichment results with NES values,
#'           statistical significance, and formatted pathway names with line
#'           breaks}
#'     \item{\code{Table_lollipopPlot}}{A tibble/data.frame containing the plot
#'           data with columns: \code{original_name} (original pathway IDs),
#'           \code{Description} (formatted pathway descriptions with line
#'           breaks),
#'           \code{NES} (Normalized Enrichment Score), \code{p.adjust} (adjusted
#'           p-values), and \code{Count} (number of genes in each pathway)}
#'   }
#' @export
#' @examples
#' data("demo_GSE181574")
#' lollipopPlot(SE_data.fgsea=GSE181574.sigfun)
lollipopPlot <- function(SE_data.fgsea, showCategory=10, breaklineN=30){
    res.gsea <- .extractDF(SE_data.fgsea, type="gseaRaw")
    if (is.numeric(showCategory)) {
        res.gsea@result <- res.gsea@result[seq_len(showCategory), ]
    }else{
        res.gsea@result <- res.gsea@result |>
            dplyr::filter(ID %in% showCategory)
    }
    lollipopPlot <- GseaVis::dotplotGsea(
        data=res.gsea, order.by='NES', add.seg=TRUE)$plot +
        ggplot2::scale_y_discrete(labels=function(x){
            .labelBreak(x, breaklineN)})
    Table_lollipopPlot <- lollipopPlot$data |>
        dplyr::mutate(Description=.labelBreak(Description, breaklineN)) |>
        dplyr::select(
            original_name=ID, Description, NES, p.adjust, Count) |>
        tibble::remove_rownames()
    return(list(
        lollipopPlot=lollipopPlot, Table_lollipopPlot=Table_lollipopPlot))
}
