#' @title emapPlot
#' @description The \code{emapPlot} function generates a network-based enrichment
#'   map visualization for Gene Set Enrichment Analysis (GSEA) results stored in
#'   a \code{SummarizedExperiment} object. Each node represents a pathway or
#'   gene set, and edges represent pairwise similarity between pathways based on
#'   shared genes. This visualization provides an overview of functional
#'   relationships among enriched terms.
#'
#' @param seDataFgsea A \code{SummarizedExperiment} object containing GSEA
#'   results and pairwise pathway similarity information, typically extracted
#'   using \code{.extractDF(type = "gseaSimilar")}.
#' @param showCategory Numeric or character. The number of top pathways to
#'   display, or a vector of specific pathway names to include. Default is
#'   \code{10}.
#' @param layout A layout function from the \pkg{igraph} package to position
#'   nodes in the network. Default is \code{igraph::layout_with_kk}.
#' @param color Character. Variable used for node color encoding. Must be one of
#'   \code{"pvalue"}, \code{"p.adjust"}, or \code{"qvalue"}. Default is
#'   \code{"p.adjust"}.
#' @param sizeCategory Numeric. Scaling factor for node sizes. Default is
#'   \code{1}.
#' @param minEdge Numeric. Minimum similarity threshold for retaining edges in
#'   the network (range: 0–1). Default is \code{0.2}.
#' @param fontSize Numeric. Font size of node labels. Default is \code{3}.
#' @param colorEdge Character. Color of connecting edges between nodes. Default
#'   is \code{"grey"}.
#' @param sizeEdge Numeric. Edge line width. Default is \code{0.4}.
#' @param nodeLabel Character. Determines which nodes are labeled in the plot.
#'   Accepts \code{"category"}, \code{"all"}, or \code{"none"}. Default is
#'   \code{"category"}.
#'
#' @details
#'   This function provides a simplified implementation of
#'   \code{enrichplot::emapplot()} with a clean SigFun-style aesthetic. Nodes
#'   are colored by the statistical significance (e.g., \code{p.adjust}) and
#'   sized according to pathway gene counts.
#'
#'   The edge structure between pathways is derived from pairwise similarity
#'   scores precomputed by \code{.extractDF(type = "gseaSimilar")}.
#'
#'   Visual encoding follows the SigFun design standard:
#'   \itemize{
#'     \item \strong{Node color:} Encodes -log10(significance)
#'     \item \strong{Node size:} Reflects pathway or gene set size
#'     \item \strong{Edge thickness:} Represents similarity between pathways
#'     \item \strong{Labels:} Optional display controlled by \code{nodeLabel}
#'   }
#'
#' @return A named list with two components:
#'   \describe{
#'     \item{\code{emapPlot}}{A \code{ggplot2} network visualization showing
#'       pathways (nodes) and their functional relationships (edges).}
#'     \item{\code{tableEmapPlot}}{A \code{data.frame} summarizing node
#'       attributes used in the plot, including:
#'       \itemize{
#'         \item \code{labelName}: Pathway or term label
#'         \item \code{size}: Node size values
#'         \item \code{color}: Color aesthetic value (when applicable)
#'       }}
#'   }
#'
#' @importFrom igraph layout_with_kk V
#' @importFrom ggtangle ggplot geom_edge %<+%
#' @importFrom ggplot2 geom_point coord_equal scale_size scale_color_gradient
#'   guides guide_colorbar guide_legend aes element_blank
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr select
#' @importFrom rlang .data
#'
#' @export
#' @examples
#' data("sig2Fun_result")
#' emapPlot(seDataFgsea = sig2Fun_result)
#'
#' @note
#' Requires \pkg{enrichplot} (>= 1.21.0) and \pkg{ggtangle} (>= 1.0.0)
#' for consistent network visualization behavior.

# Note: Requires enrichplot (>= 1.21.0) and ggtangle (>= 1.0.0)
emapPlot <- function(
        seDataFgsea, showCategory = 10,
        layout = igraph::layout_with_kk, color = "p.adjust", sizeCategory = 1,
        minEdge = .2, fontSize = 3,
        colorEdge = "grey", sizeEdge = .4,
        nodeLabel = "category"
) {

    match.arg(color, c("pvalue", "p.adjust", "qvalue"))
    x <- .extractDF(seDataFgsea, type = "gseaSimilar")

    gg <- enrichplot:::graph_from_enrichResult(
        x, showCategory = showCategory, color = color,
        min_edge = minEdge, size_edge = sizeEdge
    )

    g <- gg$graph
    size <- vapply(gg$geneSet, length, FUN.VALUE = numeric(1))
    igraph::V(g)$size <- size[match(igraph::V(g)$name, names(size))]

    p <- ggplot(g, layout = layout) +
        ggtangle::geom_edge(color = colorEdge, linewidth = sizeEdge)

    # SigFun-style aesthetic mapping
    if (color %in% names(as.data.frame(x))) {
        p <- p %<+% x[, c("Description", color)] +
            geom_point(aes(color = -log10(.data[[color]]), size = .data$size)) +
            scale_size(range = c(3, 8) * sizeCategory) +
            scale_color_gradient(
                low = "#fcae91", high = "#b3200a",
                name = bquote(-log[10](.(color))),
                labels = function(x) sprintf("%.1f", x)
            ) +
            guides(
                size = guide_legend(order = 1, override.aes = list(color = "grey60")),
                color = guide_colorbar(order = 2, reverse = TRUE, barwidth = 0.9, barheight = 4.5)
            )
    } else {
        p <- p %<+% x[, "Description", drop = FALSE] +
            geom_point(aes(size = .data$size), color = color) +
            scale_size(range = c(3, 8) * sizeCategory)
    }

    if (nodeLabel %in% c("all", "category")) {
        p <- p + ggrepel::geom_text_repel(
            aes(label = .data$label), bg.color = "white", bg.r = .1, size = fontSize
        )
    }

    p <- p + coord_equal()
    tableEmapPlot <- p$data |>
        dplyr::select(labelName = label, size, color)
    return(list(emapPlot = p, tableEmapPlot = tableEmapPlot))
}
