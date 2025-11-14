#' @title treePlot
#' @description The \code{treePlot} function creates a hierarchical tree plot
#'   to visualize gene set enrichment analysis (GSEA) results. The tree plot
#'   groups similar pathways together based on gene overlap and displays them
#'   as a dendrogram, providing an intuitive way to explore functional clusters
#'   and pathway relationships. Color and size mappings represent enrichment
#'   significance and gene counts, respectively.
#'
#' @param seDataFgsea A \code{SummarizedExperiment} object containing GSEA
#'   results in its metadata. The object must include:
#'   \itemize{
#'     \item \code{metadata(seDataFgsea)$gseaResult}: GSEA result object with
#'           pathway information, enrichment statistics, and core enrichment genes.
#'   }
#' @param showCategory Integer or character. Specifies the number of top
#'   pathways to display, or a vector of pathway IDs/names to include.
#'   Must include at least 5 pathways. Default is \code{10}.
#' @param dotColor Character. Column name in the GSEA results used to map color.
#'   Options: \code{'pvalue'}, \code{'p.adjust'}, \code{'qvalue'}, or \code{'NES'}.
#'   Default is \code{'pvalue'}.
#' @param hclustfun Character. Hierarchical clustering method used to group
#'   pathways. Options include \code{'ward.D'}, \code{'ward.D2'},
#'   \code{'single'}, \code{'complete'}, \code{'average'}, \code{'mcquitty'},
#'   \code{'median'}, or \code{'centroid'}. Default is \code{'ward.D'}.
#' @param labelFormat Optional. Numeric or function controlling text wrapping
#'   or formatting for clade labels.
#' @param labelFormatTiplab Optional. Numeric or function controlling text
#'   formatting for terminal node (tip) labels.
#' @param cexCategory Numeric. Scaling factor for node size and label spacing.
#'   Default is \code{1}.
#' @param leafFontSize Numeric. Font size for terminal node labels.
#'   Default is \code{3}.
#' @param cladeFontSize Numeric. Font size for clade (internal node) labels.
#'   Default is \code{2.5}.
#' @param hilightParams List. Controls highlighting and alignment behavior for
#'   clades. Default is \code{list(hilight = TRUE, align = "both")}.
#'   \itemize{
#'     \item \code{hilight}: Whether to shade cluster regions.
#'     \item \code{align}: Text alignment mode for cluster labels.
#'   }
#' @param offsetParams List. Controls tree and label positioning parameters.
#'   Default is \code{list(barTree = rel(1.3), tiplab = rel(1.5),
#'   extend = 0.3, hexpand = 0.1)}.
#'   \itemize{
#'     \item \code{barTree}: Horizontal spacing between tree branches and text.
#'     \item \code{tiplab}: Distance between tips and labels.
#'     \item \code{extend}: Extension of the tree layout for clarity.
#'     \item \code{hexpand}: Expansion ratio applied to x-axis layout.
#'   }
#' @param clusterParams List. Parameters for hierarchical clustering and labeling.
#'   Default is \code{list(method = "ward.D", n = 5, color = NULL,
#'   labelWordsN = 0, labelFormat = 30)}.
#'   \itemize{
#'     \item \code{method}: Clustering algorithm.
#'     \item \code{n}: Number of clusters to divide the tree into.
#'     \item \code{color}: Optional vector for custom cluster colors.
#'     \item \code{labelWordsN}: Maximum number of words shown in clade labels.
#'     \item \code{labelFormat}: Text wrapping width for clade labels.
#'   }
#'
#' @details
#'   The \code{treePlot} function organizes enriched pathways into functional
#'   clusters based on term similarity (gene overlap) and visualizes them using
#'   a dendrogram structure. Users can customize appearance through nested
#'   parameter lists (\code{hilightParams}, \code{offsetParams},
#'   \code{clusterParams}) to fine-tune alignment, spacing, and label style.
#'
#'   Internal hierarchical clustering is computed via \code{stats::hclust},
#'   and visualization is rendered using \code{ggtree}. Terminal nodes (tips)
#'   correspond to enriched pathways, with color intensity representing
#'   enrichment significance and point size corresponding to the number of
#'   core genes per pathway.
#'
#' @return A list containing two elements:
#'   \describe{
#'     \item{\code{treePlot}}{A \code{ggtree} object displaying the hierarchical
#'       clustering of pathways with color-coded enrichment metrics and cluster
#'       highlights.}
#'     \item{\code{tableTreePlot}}{A \code{data.frame} containing pathway-level
#'       data for terminal nodes, including label names, selected color variable
#'       (e.g., \code{pvalue}, \code{p.adjust}, or \code{NES}), and gene counts.}
#'   }
#'
#' @importFrom ggtree ggtree MRCA groupOTU hexpand
#' @importFrom ggplot2 scale_size_continuous guide_legend guide_colorbar
#' @importFrom stats hclust cutree as.dist
#' @importFrom dplyr select filter
#' @importFrom cli cli_abort
#' @export
#'
#' @examples
#' data("sig2Fun_result")
#' treePlot(seDataFgsea = sig2Fun_result,
#'           showCategory = 10,
#'           dotColor = "pvalue",
#'           hclustfun = "ward.D")
treePlot <- function(seDataFgsea, showCategory = 10,
                     dotColor = 'pvalue', hclustfun = 'ward.D',
                     labelFormat = NULL, labelFormatTiplab = NULL,
                     cexCategory = 1, leafFontSize = 3, cladeFontSize = 2.5,
                     hilightParams = list(hilight = TRUE, align = "both"),
                     offsetParams = list(barTree = rel(1.3), tiplab = rel(1.5), extend = 0.3, hexpand = .1),
                     clusterParams = list(method = "ward.D", n = 5, color = NULL,
                                          labelWordsN = 0, labelFormat = 30)) {

    match.arg(dotColor, c('pvalue', 'p.adjust', 'qvalue', 'NES'))
    x <- .extractDF(seDataFgsea, type = "gseaSimilar")

    if (is.numeric(showCategory)) {
        if (showCategory < 5) cli::cli_abort("{.arg showCategory} must be larger than {.val 5}.")
    } else {
        if (sum(showCategory %in% x@result$ID) < 5)
            cli::cli_abort("{.arg showCategory} must provide at least {.val 5} pathway names.")
    }

    paramsDf <- as.data.frame(rbind(
        c("hilight", "hilightParams", "hilight"),
        c("align", "hilightParams", "align"),
        c("offset", "offsetParams", "barTree"),
        c("offsetTiplab", "offsetParams", "tiplab"),
        c("extend", "offsetParams", "extend"),
        c("hexpand", "offsetParams", "hexpand"),
        c("hclustMethod", "clusterParams", "method"),
        c("nCluster", "clusterParams", "n"),
        c("groupColor", "clusterParams", "color"),
        c("nWords", "clusterParams", "labelWordsN"),
        c("labelFormatCladelab", "clusterParams", "labelFormat"))
    )
    colnames(paramsDf) <- c("original", "listname", "present")
    rownames(paramsDf) <- paramsDf$original

    # Default parameter sets
    defaultHilightParams <- list(hilight = TRUE, align = "both")
    defaultOffsetParams  <- list(barTree = rel(1.5), tiplab = rel(1.5), extend = 0.3, hexpand = .1)
    defaultClusterParams <- list(method = hclustfun, n = 5, color = NULL, labelWordsN = 0, labelFormat = 30)

    hilightParams <- modifyList(defaultHilightParams, hilightParams)
    offsetParams  <- modifyList(defaultOffsetParams, offsetParams)
    clusterParams <- modifyList(defaultClusterParams, clusterParams)
    # Extract essential parameters from clusterParams before building paramsList
    nWords <- clusterParams$labelWordsN
    nCluster <- clusterParams$n
    groupColor <- clusterParams$color
    labelFormatCladelab <- clusterParams$labelFormat
    hclustMethod <- clusterParams$method
    paramsList <- list(
        x = x,
        showCategory = showCategory,
        dotColor = dotColor,
        nWords = nWords,
        nCluster = nCluster,
        cexCategory = cexCategory,
        labelFormat = labelFormat,
        labelFormatCladelab = labelFormatCladelab,
        labelFormatTiplab = labelFormatTiplab,
        leafFontSize = leafFontSize,
        cladeFontSize = cladeFontSize,
        offset = offsetParams$barTree,
        offsetTiplab = offsetParams$tiplab,
        hclustMethod = hclustMethod,
        groupColor = groupColor,
        extend = offsetParams$extend,
        hilight = hilightParams$hilight,
        hexpand = offsetParams$hexpand,
        align = hilightParams$align,
        hilightParams = hilightParams,
        offsetParams = offsetParams,
        clusterParams = clusterParams
    )

    args <- as.list(match.call())
    removedParams <- intersect(paramsDf$original, names(args))
    if (length(removedParams) > 0) {
        for (i in removedParams) {
            paramsList[[paramsDf[i, 2]]][[paramsDf[i, 3]]] <- get(i)
            warn <- getParamChangeMessage(i, paramsDf)
            warning(warn)
        }
    }

    hilightParams <- paramsList[["hilightParams"]]
    offsetParams  <- paramsList[["offsetParams"]]
    clusterParams <- paramsList[["clusterParams"]]

    hilight <- hilightParams$hilight
    align <- hilightParams$align
    offset <- offsetParams$barTree
    offsetTiplab <- offsetParams$tiplab
    extend <- offsetParams$extend
    hexpand <- offsetParams$hexpand
    hclustMethod <- clusterParams$method
    nCluster <- clusterParams$n
    groupColor <- clusterParams$color
    nWords <- clusterParams$labelWordsN
    labelFormatCladelab <- clusterParams$labelFormat

    group <- p.adjust <- count <- NULL
    if (!is.null(labelFormat)) labelFormatCladelab <- labelFormat

    if (inherits(x, "gseaResult")) {
        x@result$Count <- x$core_enrichment %>%
            strsplit(split = "/") %>%
            vapply(length, FUN.VALUE = 1)
    }

    n <- enrichplot:::update_n(x, showCategory)
    keep <- if (is.numeric(n)) seq_len(n) else match(n, rownames(x@termsim))
    if (length(keep) == 0) stop("no enriched term found...")

    termSim2 <- enrichplot:::fill_termsim(x, keep)
    hc <- stats::hclust(stats::as.dist(1 - termSim2), method = hclustMethod)
    clus <- stats::cutree(hc, nCluster)
    d <- data.frame(
        label = names(clus),
        color = as.numeric(x[keep, as.character(dotColor)]),
        count = x$Count[keep]
    )

    # plotting
    treePlot <- suppressWarnings(.groupTree(
        hc = hc,
        clus = clus,
        d = d,
        offsetTiplab = offsetTiplab,
        nWords = nWords,
        labelFormatCladelab = labelFormatCladelab,
        labelFormatTiplab = labelFormatTiplab,
        offset = offset,
        leafFontSize = leafFontSize,
        cladeFontSize = cladeFontSize,
        groupColor = groupColor,
        extend = extend,
        hilight = hilight,
        cexCategory = cexCategory,
        align = align,
        dotColor = dotColor
    ))

    treePlot <- treePlot +
        ggplot2::scale_size_continuous(name = "number of genes", range = c(2, 5) * cexCategory) +
        ggtree::hexpand(ratio = hexpand) +
        ggplot2::guides(
            size = ggplot2::guide_legend(order = 1, override.aes = list(color = "grey60")),
            color = ggplot2::guide_colorbar(order = 2, barwidth = 0.9, barheight = 4.5)
        )

    tableTreePlot <- treePlot$data %>%
        dplyr::filter(isTip == TRUE) %>%
        dplyr::select(labelName = label, !!dotColor := color, count)
    return(list(treePlot = treePlot, tableTreePlot = tableTreePlot))
}
