#' @title emapPlot
#' @description The \code{emapPlot} function creates an enrichment map plot to
#' visualize gene set enrichment analysis (GSEA) results. This network-based
#' visualization displays pathways as nodes and connects similar pathways with
#' edges, where edge thickness represents the degree of similarity between
#' pathways. The enrichment map provides an intuitive overview of the
#' functional landscape by clustering related biological processes and
#' highlighting pathway relationships through shared genes.
#' @param SE_data.fgsea A \code{SummarizedExperiment} object containing GSEA
#' results with precomputed pairwise term similarities. The object must contain
#' processed GSEA data that can be extracted using the internal helper function
#' \code{.extractDF} with \code{type="gseaSim"}. This typically includes:
#'   \itemize{
#'     \item GSEA results with readable gene symbols
#'     \item Pairwise pathway similarity calculations
#'     \item Enrichment statistics and significance values
#'   }
#' @param showCategory An integer specifying the maximum number of pathways to
#' display in the enrichment map. Pathways are typically ranked by significance
#' or enrichment score, and the top pathways are selected for visualization.
#' Default is \code{30}.
#' @return A list containing two elements:
#'   \describe{
#'     \item{\code{emapPlot}}{A ggplot2 object containing the enrichment map
#'           visualization showing pathway networks with nodes representing
#'           pathways and edges representing functional similarities between
#'           pathways}
#'     \item{\code{Table_emapPlot}}{A data.frame containing the plot data
#'           with columns: \code{label_name} (pathway names/labels),
#'           \code{size} (node size values, typically representing gene count
#'           or significance), and \code{p.adjust} (adjusted p-values for
#'           statistical significance)}
#'   }
#' @export
#' @examples
#' data("demo_GSE181574")
#' emapPlot(SE_data.fgsea=GSE181574.sigfun)
emapPlot <- function(SE_data.fgsea, showCategory=10){
    gseaSim <- .extractDF(SE_data.fgsea, type="gseaSimilar")
    emapPlot <- enrichplot::emapplot(gseaSim, showCategory=showCategory)
    Table_emapPlot <- emapPlot$data |>
        dplyr::select(label_name=label, size, p.adjust)
    return(list(emapPlot=emapPlot, Table_emapPlot=Table_emapPlot))
}
