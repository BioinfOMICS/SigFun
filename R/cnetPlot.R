#' @title cnetPlot
#'
#' @description The \code{cnetPlot} function generates a highly customizable
#' concept network (CNET) visualization for gene set enrichment analysis (GSEA)
#' results. It displays the relationships between enriched pathways and their
#' associated genes as a network graph, where nodes represent pathways or genes
#' and edges represent their associations. Users can control node and edge
#' aesthetics, highlight specific pathways, and optionally map gene-level
#' correlation values to node colors.
#'
#' @param seDataFgsea A \code{SummarizedExperiment} object containing GSEA
#'   results. The object must include:
#'   \itemize{
#'     \item \code{metadata(seDataFgsea)$gseaResult}: A GSEA result object with
#'           pathway information, gene sets, and enrichment statistics.
#'     \item (Optional) \code{metadata(seDataFgsea)$cor.df}: A data frame of gene-level
#'           correlation values, required when using the \code{color} parameter.
#'   }
#' @param showCategory Either a numeric value specifying the number of top
#'   pathways to display, or a character vector specifying specific pathway
#'   names. Default is \code{5}.
#' @param nodeLabel Character string specifying which type of nodes should be
#'   labeled. Must be one of \code{c("category", "all", "none", "item", "gene",
#'   "exclusive", "share")}. Default is \code{"category"}.
#' @param fontSize Numeric. Base size for text labels (applied as a scaling factor).
#'   Default is \code{1}.
#' @param color Character or \code{NULL}. The column name from
#'   \code{metadata(seDataFgsea)$cor.df} to be used for coloring gene nodes
#'   (e.g., correlation coefficients or log2 fold changes). When \code{NULL},
#'   nodes are colored using the default item color. Default is \code{NULL}.
#' @param circular Logical. Whether to arrange nodes in a circular layout.
#'   Default is \code{FALSE}.
#' @param breaklineN Integer. Maximum number of characters before inserting line
#'   breaks in pathway names for improved readability. Default is \code{30}.
#' @param layout Layout function for node positioning. Default is
#'   \code{igraph::layout_nicely}. Can be replaced with any \pkg{igraph} layout
#'   function.
#' @param colorCategory Character. Color for category (pathway) nodes.
#'   Default is \code{"#E5C494"}.
#' @param sizeCategory Numeric. Scaling factor for category node sizes.
#'   Default is \code{0.9}.
#' @param colorItem Character. Color for item (gene) nodes when
#'   \code{color = NULL}. Default is \code{"#B3B3B3"}.
#' @param sizeItem Numeric. Scaling factor for item (gene) node sizes.
#'   Default is \code{0.8}.
#' @param colorEdge Character. Edge color or \code{"category"} to color edges by
#'   their associated pathway. Default is \code{"grey80"}.
#' @param sizeEdge Numeric. Edge line width. Default is \code{0.4}.
#' @param hilight Character vector or \code{"none"}. Pathways to highlight
#'   (opacity = 1). Non-highlighted nodes will have reduced opacity controlled
#'   by \code{hilightAlpha}. Default is \code{"none"}.
#' @param hilightAlpha Numeric. Transparency level (0-1) for non-highlighted
#'   nodes. Default is \code{0.3}.
#' @param seed Integer or \code{NULL}. Random seed for reproducible layout
#'   generation. Default is \code{NULL}.
#'
#' @details
#' This function builds upon \pkg{ggtangle} and \pkg{igraph} for flexible
#' concept network visualization. Users can:
#' \itemize{
#'   \item Control node/edge colors, sizes, and layout.
#'   \item Overlay continuous gene metrics (e.g., correlation, fold change).
#'   \item Highlight specific pathways using \code{hilight}.
#'   \item Automatically wrap long pathway names using \code{breaklineN}.
#' }
#'
#' When \code{color} is provided, gene nodes are colored according to their
#' corresponding values using a diverging or sequential color scale. When
#' \code{color = NULL}, all item nodes are rendered in a uniform grey tone.
#'
#' @return A named list containing three components:
#' \itemize{
#'   \item \code{cnetPlot}: A \code{ggplot2} object representing the concept
#'         network visualization.
#'   \item \code{tableEdgeCnetPlot}: A \code{data.frame} of network edges with
#'         columns:
#'         \itemize{
#'           \item \code{from}: Pathway (category) identifier.
#'           \item \code{to}: Gene symbol(s) associated with the pathway.
#'         }
#'   \item \code{tableNodeCnetPlot}: A \code{data.frame} of node attributes with
#'         columns:
#'         \itemize{
#'           \item \code{id}: Node identifier (pathway or gene name).
#'           \item \code{size}: Node size used in visualization.
#'           \item (Optional) Column named by \code{color} parameter, containing
#'                 the numeric values mapped to gene color scale.
#'         }
#' }
#'
#' @importFrom ggtangle geom_cnet_label td_filter
#' @importFrom igraph layout_nicely V graph_from_data_frame degree
#' @importFrom igraph as_edgelist edge_attr edge_attr_names
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient scale_color_gradient2
#' @importFrom dplyr select mutate pull
#' @importFrom tidyr separate_rows
#' @importFrom cli cli_abort cli_alert_warning
#'
#' @export
#' @examples
#' data("sig2Fun_result")
#' cnetPlot(seDataFgsea = sig2Fun_result,
#'          color = "cor")
#'
#' @note
#' Requires \pkg{ggtangle} (>= 1.0.0), \pkg{enrichplot} (>= 1.21.0),
#' and \pkg{ggplot2} (>= 3.5.0) for full compatibility.

# Note: Requires ggtangle (>= 1.0.0), enrichplot (>= 1.21.0), and ggplot2 (>= 3.5.0)
cnetPlot <- function(
        seDataFgsea, showCategory = 5, nodeLabel = "category",
        fontSize = 1, color = NULL, circular = FALSE, breaklineN = 30,
        layout = igraph::layout_nicely,
        colorCategory = "#E5C494", sizeCategory = 0.9,
        colorItem = "#B3B3B3",  sizeItem  = 0.8,
        colorEdge = "grey80",   sizeEdge  = 0.4,
        hilight = "none",       hilightAlpha = 0.3,
        seed = NULL
) {

    gseaObj <- .extractDF(seDataFgsea, type = "gseaReadable")

    if (!is.null(color)) {
        corDf <- .extractDF(seDataFgsea, type = "corCoef")
        if (gseaObj@readable && gseaObj@keytype != "SYMBOL") {
            corDf <- data.frame(gene = names(gseaObj@gene2Symbol),
                                geneSymbol = gseaObj@gene2Symbol) |>
                merge(corDf) |>
                dplyr::select(-gene) |>
                dplyr::rename(gene = geneSymbol)
        }
        if (!(color %in% colnames(corDf))) {
            cli::cli_abort(
                "The input {.arg color} must be a {.cls character} and must
                be included in the column names of
                {.field seDataFgsea@metadata$cor.df}."
            )
        }
        foldChange <- corDf |>
            dplyr::pull(.data[[color]]) |>
            stats::setNames(corDf$gene)
    } else {
        foldChange <- NULL
    }

    x <- DOSE::geneInCategory(gseaObj)[gseaObj@result$ID]
    origNames <- names(x)  # store original (unwrapped) names

    if (is.numeric(showCategory)) {
        x <- x[seq_len(showCategory)]
    } else {
        matchedIdx <- match(showCategory, names(x))
        matchedIdx <- matchedIdx[!is.na(matchedIdx)]
        if (length(matchedIdx) == 0) {
            cli::cli_abort("No matching pathways found in {.field seDataFgsea@result}.")
        }
        x <- x[matchedIdx]
    }

    # for hilight matching
    wrappedNames <- .labelBreak(names(x), breaklineN)
    attr(x, "mapWrappedToOrig") <- stats::setNames(names(x), wrappedNames)
    names(x) <- wrappedNames

    if (is.numeric(showCategory)) {
        x <- .subsetCnetList(x, showCategory)
    } else {
        # Already handled by the match() function above
        x <- x
    }
    cnt <- stats::setNames(sapply(x, length), names(x))

    if (length(nodeLabel) > 1) {
        if (getOption("cnetplot_subset", default = FALSE)) {
            x <- .subsetCnetListItem(x, nodeLabel)
            nodeLabel <- "all"
        }
    } else if (!nodeLabel %in% c("category","all","none","item","gene","exclusive","share")) {
        if (!grepl("[><=]", nodeLabel)) {
            stop("wrong parameter for 'nodeLabel'")
        } else if (is.null(foldChange)) {
            stop("'foldChange' should not be NULL with the 'nodeLabel' setting")
        }
    }
    if (length(nodeLabel) == 1 && nodeLabel == "gene") nodeLabel <- "item"

    g <- .list2graph(x)

    igraph::V(g)$`.hilight` <- 1
    # Use mapWrappedToOrig to ensure correct matching after .labelBreak()
    if (all(hilight != "none")) {
        mapWrappedToOrig <- attr(x, "mapWrappedToOrig")
        if (!is.null(mapWrappedToOrig)) {
            matchWrapped <- names(mapWrappedToOrig)[mapWrappedToOrig %in% hilight |
                                                        names(mapWrappedToOrig) %in% hilight]
            if (length(matchWrapped) > 0) {
                y <- .subsetCnetList(x, matchWrapped)
                igraph::V(g)$`.hilight` <- hilightAlpha
                igraph::V(g)$`.hilight`[igraph::V(g)$name %in% names(y)] <- 1
                igraph::V(g)$`.hilight`[igraph::V(g)$name %in% unlist(y)] <- 1
            } else {
                cli::cli_alert_warning("No matching hilight pathways found - all nodes remain default opacity.")
            }
        }
    }

    if (!is.null(foldChange)) {
        igraph::V(g)$foldChange <- foldChange[igraph::V(g)$name]
        fcMapping <- ggplot2::aes(color = .data$foldChange, alpha = I(.data$.hilight))
    } else {
        fcMapping <- ggplot2::aes(color = I(colorItem), alpha = I(.data$.hilight))
    }

    #set.seed
    if (!is.null(seed)) set.seed(seed)

    p <- ggtangle::ggplot(g, layout = layout)

    ## restore original category size
    if (length(nodeLabel) > 1 &&
        getOption("cnetplot_subset", default = FALSE)) {
        p$data$size[match(names(cnt), p$data$label)] <- cnt
    }

    if (colorEdge == "category") {
        ed <- .getEdgeData(g)
        names(ed)[1] <- "category"
        p <- p + ggtangle::geom_edge(ggplot2::aes(color = .data$category), data = ed, linewidth = sizeEdge) +
            ggnewscale::new_scale_color()
    } else {
        p <- p + ggtangle::geom_edge(color = colorEdge, linewidth = sizeEdge)
    }

    p <- p + ggplot2::geom_point(ggplot2::aes(size = .data$size, alpha = I(.data$.hilight)),
                                 data = ggtangle::td_filter(.data$.isCategory),
                                 color = colorCategory) +
        ggplot2::geom_point(fcMapping,
                            data = ggtangle::td_filter(!.data$.isCategory), size = 3 * sizeItem) +
        ggplot2::scale_size(range = c(3, 8) * sizeCategory, breaks = pretty(cnt, n = min(4, diff(range(cnt)))))

    # fontSize parameter to control the text size
    if (length(nodeLabel) > 1 || nodeLabel != "none") {
        if (length(nodeLabel) > 1 ||
            nodeLabel %in% c("exclusive","share")) {
            p <- p + ggtangle::geom_cnet_label(node_label = "category", size = fontSize * 3)
        }
        p <- p + ggtangle::geom_cnet_label(node_label = nodeLabel, size = fontSize * 3)
    }

    # SigFun-style aesthetic mapping
    if (!is.null(foldChange)) {
        vals <- igraph::V(g)$foldChange
        rng <- range(vals, na.rm = TRUE)
        if (rng[1] >= 0) {
            p <- p + ggplot2::scale_color_gradient(
                name = color, low = "#fcae91", high = "#b3200a",
                guide = guide_colorbar(barwidth = 0.9, barheight = 4.5)
            )
        } else if (rng[2] <= 0) {
            p <- p + ggplot2::scale_color_gradient(
                name = color, low = "#08306b", high = "#c6dbef",
                guide = guide_colorbar(barwidth = 0.9, barheight = 4.5)
            )
        } else {
            p <- p + ggplot2::scale_color_gradient2(
                name = color, low = "#08306b", mid = "white", high = "#b3200a", midpoint = 0,
                guide = guide_colorbar(barwidth = 0.9, barheight = 4.5)
            )
        }
    }
    p <- p +
        ggplot2::theme(plot.margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))

    tableEdgeCnetPlot <- gseaObj@result |>
        dplyr::select(from = ID, to = core_enrichment) |>
        tidyr::separate_rows(to, sep = "/")

    tableNodeCnetPlot <- data.frame(
        id   = igraph::V(g)$name,
        size = ifelse(igraph::V(g)$name %in% names(cnt), cnt[igraph::V(g)$name], 1),
        stringsAsFactors = FALSE
    )
    if (!is.null(foldChange)) {tableNodeCnetPlot[[color]] <- foldChange[tableNodeCnetPlot$id]}

    list(
        cnetPlot = p,
        tableEdgeCnetPlot = tableEdgeCnetPlot,
        tableNodeCnetPlot = tableNodeCnetPlot
    )
}
