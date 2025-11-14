#' @title gseaPlot
#' @description The \code{gseaPlot} function generates a publication-quality
#'   visualization of Gene Set Enrichment Analysis (GSEA) results from a
#'   \code{SummarizedExperiment} object. It wraps and extends
#'   \code{enrichplot::gseaplot2()} with additional control over layout,
#'   color palette, font size, subplot arrangement, and significance table
#'   display options.
#'
#' @param seDataFgsea A \code{SummarizedExperiment} object containing GSEA
#'   results that can be processed by the internal \code{.extractDF()} function
#'   with \code{type = "gseaRaw"}.
#' @param showCategory Numeric or character. The number of top-ranked gene sets
#'   to visualize, or specific pathway names to plot. Default is \code{5}.
#' @param title Character. Optional plot title displayed above the enrichment
#'   score curve. Default is an empty string (\code{""}).
#' @param colorPalette Character vector. Custom color palette for pathways.
#'   If \code{NULL}, the default \pkg{viridisLite} palette ("D") is used.
#'   Default is \code{NULL}.
#' @param fontSize Numeric. Base font size for plot text. Default is \code{10}.
#' @param relHeights Numeric vector. Relative heights for the main enrichment
#'   score plot, gene tick mark panel, and ranked metric plot.
#'   Default is \code{c(1.5, 0.5, 1)}.
#' @param subPlots Integer vector. Specifies which subplot(s) to display:
#'   \code{1} = running score, \code{2} = gene ticks, \code{3} = ranked list.
#'   Default is \code{1:3}.
#' @param pvalueTable Logical. Whether to include a p-value summary table in the
#'   enrichment plot. Default is \code{FALSE}.
#' @param pvalueTableColumns Character vector specifying the statistical
#'   columns (e.g., \code{"pvalue"}, \code{"p.adjust"}) to include in the table.
#'   Default is \code{c("pvalue", "p.adjust")}.
#' @param ESgeom Character. Type of geometry for drawing the enrichment score
#'   curve. Must be one of \code{"line"} or \code{"dot"}. Default is
#'   \code{"line"}.
#'
#' @details
#'   This function produces a comprehensive visualization combining the
#'   running enrichment score, gene positions, and ranked list metrics in a
#'   unified layout. It improves upon \code{enrichplot::gseaplot2()} by
#'   introducing:
#'   \itemize{
#'     \item Customizable subplot arrangement (\code{subPlots})
#'     \item Optional significance table overlay (\code{pvalueTable = TRUE})
#'     \item Adjustable color palettes via \code{colorPalette}
#'     \item Flexible layout scaling using \code{relHeights}
#'     \item Support for multiple gene sets or a single enriched term
#'   }
#'
#' @return A single \code{ggplot2} object if one subplot is selected, or a
#'   combined \code{aplot::gglist} object when multiple subplots are rendered.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_linerange geom_rect
#'   geom_segment ggtitle theme_classic theme element_line element_blank
#'   element_text element_rect scale_color_viridis_d scale_color_manual
#'   scale_x_continuous scale_y_continuous xlab ylab labs annotation_custom
#'   margin
#' @importFrom viridisLite viridis
#' @importFrom aplot gglist
#' @importFrom cli cli_alert_warning
#' @importFrom stats quantile
#' @importFrom rlang .data
#'
#' @export
#' @examples
#' data("sig2Fun_result")
#' gseaPlot(seDataFgsea = sig2Fun_result)
#'
#' @note
#' Requires \pkg{ggplot2} (>= 3.5.0) for proper support of
#' \code{legend.position = "inside"} positioning.

# Note: Requires ggplot2 (>= 3.5.0) for `legend.position = "inside"` support.
gseaPlot <- function(
        seDataFgsea, showCategory = 5,
        title = "", colorPalette = NULL, fontSize = 10,
        relHeights = c(1.5, 0.5, 1), subPlots = 1:3,
        pvalueTable = FALSE, pvalueTableColumns = c("pvalue", "p.adjust"),
        ESgeom = "line"
) {

    x <- .extractDF(seDataFgsea, type = "gseaRaw")
    if (is.numeric(showCategory)) {
        geneSetID <- seq_len(showCategory)
    } else {
        geneSetID <- showCategory
    }

    ESgeom <- match.arg(ESgeom, c("line", "dot"))
    geneList <- position <- NULL
    gsdata <- enrichplot:::get_gsdata(x, geneSetID)

    p <- ggplot2::ggplot(gsdata, ggplot2::aes(x = x)) + ggplot2::xlab(NULL) +
        ggplot2::theme_classic(fontSize) +
        ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey95"),
                       panel.grid.minor = ggplot2::element_line(colour = "grey95"),
                       panel.grid.major.y = ggplot2::element_blank(),
                       panel.grid.minor.y = ggplot2::element_blank()) +
        ggplot2::scale_x_continuous(expand = c(0,0))

    if (ESgeom == "line") {
        esLayer <- ggplot2::geom_line(ggplot2::aes(y = runningScore, color = Description), linewidth = 0.8)
    } else {
        esLayer <- ggplot2::geom_point(ggplot2::aes(y = runningScore, color = Description),
                                       size = 1, data = subset(gsdata, position == 1))
    }

    p.res <- p + esLayer +
        theme(legend.position="inside",
              legend.position.inside = c(.8, .8), legend.title = element_blank(),
              legend.background = element_rect(fill = "transparent"))

    p.res <- p.res + ylab("Running Enrichment Score") +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.line.x=element_blank(),
              plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))

    i <- 0
    for (term in unique(gsdata$Description)) {
        idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
        gsdata[idx, "ymin"] <- i
        gsdata[idx, "ymax"] <- i + 1
        i <- i + 1
    }
    p2 <- ggplot2::ggplot(gsdata, ggplot2::aes(x = x)) +
        ggplot2::geom_linerange(ggplot2::aes(ymin = ymin, ymax = ymax, color = Description)) +
        ggplot2::xlab(NULL) + ggplot2::ylab(NULL) + ggplot2::theme_classic(fontSize) +
        ggplot2::theme(legend.position = "none",
                       plot.margin = ggplot2::margin(t = -0.1, b = 0, unit = "cm"),
                       axis.ticks = ggplot2::element_blank(),
                       axis.text = ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank()) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0))

    if (length(geneSetID) == 1) {
        v <- seq(1, sum(gsdata$position), length.out = 9)
        inv <- findInterval(rev(seq_along(gsdata$position)), v)
        if (min(inv) == 0) inv <- inv + 1

        # viridis option "D" color palette
        col <- viridisLite::viridis(9, option = "D")

        ymin <- min(p2$data$ymin)
        yy <- max(p2$data$ymax - p2$data$ymin) * .3
        xmin <- which(!duplicated(inv))
        xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
        d <- data.frame(ymin = ymin, ymax = yy,
                        xmin = xmin,
                        xmax = xmax,
                        col = col[unique(inv)])
        # Updated for modern ggplot2 syntax
        p2 <- p2 + geom_rect(
            aes(
                xmin = .data$xmin,
                xmax = .data$xmax,
                ymin = .data$ymin,
                ymax = .data$ymax,
                fill = I(.data$col)
            ),
            data = d,
            alpha = 0.9,
            inherit.aes = FALSE
        )
    }

    df2 <- p$data
    df2$y <- p$data$geneList[df2$x]
    p.pos <- p + geom_segment(
        data = df2, aes(x = .data$x, xend = .data$x, y = .data$y, yend = 0), color = "grey")
    p.pos <- p.pos + ylab("Ranked List Metric") +
        xlab("Rank in Ordered Dataset") +
        theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))

    # plot title with bold font and increased size if title is provided
    if (!is.null(title) && !is.na(title) && title != "")
        p.res <- p.res +
        ggplot2::ggtitle(title) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = fontSize + 4, face = "bold")
        )

    # colorPalette parameter setting
    if (is.null(colorPalette)) {
        p.res <- p.res + ggplot2::scale_color_viridis_d(option = "D")
        p2    <- p2    + ggplot2::scale_color_viridis_d(option = "D")
    } else {
        n_path <- length(unique(gsdata$Description))
        n_col  <- length(colorPalette)
        if (n_col != n_path) {
            cli::cli_alert_warning(
                paste0("`colorPalette` length does not match number of pathways. Using default viridis palette instead.")
            )
            p.res <- p.res + ggplot2::scale_color_viridis_d(option = "D")
            p2    <- p2    + ggplot2::scale_color_viridis_d(option = "D")
        } else {
            p.res <- p.res + ggplot2::scale_color_manual(values = colorPalette)
            p2    <- p2    + ggplot2::scale_color_manual(values = colorPalette)
            if (n_path == 1) {
                p.res <- p.res + ggplot2::theme(legend.position = "none")
            }
        }
    }

    # The pvalueTableRownames parameter is no longer used and has been removed
    if (pvalueTable) {
        pd <- x[geneSetID, pvalueTableColumns]
        for (i in seq_len(ncol(pd))) {
            pd[, i] <- format(pd[, i], digits = 4)
        }
        # Replace deprecated `rows` argument with default behavior
        tp <- enrichplot:::tableGrob2(d = pd, p = p.res)

        # Reduce font size of all text grobs in the table for better readability
        text_idx <- which(sapply(tp$grobs, inherits, "text"))
        for (i in text_idx) {
            tp$grobs[[i]]$gp$fontsize <- max(fontSize - 2, 6)
        }


        p.res <- p.res + theme(legend.position = "none") +
            annotation_custom(tp,
                              xmin = quantile(p.res$data$x, .5),
                              xmax = quantile(p.res$data$x, .95),
                              ymin = quantile(p.res$data$runningScore, .7),
                              ymax = quantile(p.res$data$runningScore, .9))
    }

    plotlist <- list(p.res, p2, p.pos)[subPlots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] +
        theme(axis.line.x = element_line(),
              axis.ticks.x=element_line(),
              axis.text.x = element_text())

    if (length(subPlots) == 1)
        return(plotlist[[1]] + theme(plot.margin=margin(t=.2, r = .2, b=.2,
                                                        l=.2, unit="cm")))


    if (length(relHeights) > length(subPlots))
        relHeights <- relHeights[subPlots]

    aplot::gglist(gglist = plotlist, ncol=1, heights=relHeights)

}
