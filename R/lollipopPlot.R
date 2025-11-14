#' @title lollipopPlot
#' @description The \code{lollipopPlot} function generates a lollipop-style
#'   visualization of Gene Set Enrichment Analysis (GSEA) results. Each
#'   horizontal line represents a pathway, and the dot at its end marks the
#'   enrichment statistic (e.g., NES or GeneRatio). Dot size and transparency
#'   encode statistical significance, while color indicates enrichment
#'   direction. The connecting line segments can be toggled on or off.
#'
#' @param seDataFgsea A \code{SummarizedExperiment} object containing GSEA
#'   results in its metadata.
#' @param showCategory Numeric or character. Number of top pathways to display
#'   or specific pathway names to include. When numeric:
#'   \itemize{
#'     \item If \code{x = "NES"}: Displays top N positively enriched and top N
#'       negatively enriched pathways (up to 2N total pathways)
#'     \item If \code{x} is \code{"Count"} or \code{"GeneRatio"}: Displays top N
#'       pathways sorted by the specified metric
#'   }
#'   When character: Displays only the specified pathway IDs/names.
#'   Default is \code{10}.
#' @param breaklineN Integer. Maximum number of characters before wrapping
#'   pathway names. Default is \code{30}.
#' @param x Character. Variable for the x-axis. Must be one of \code{"NES"},
#'   \code{"Count"}, or \code{"GeneRatio"}. Default is \code{"NES"}.
#' @param alpha Character. Variable for transparency encoding. Must be one of
#'   \code{"pvalue"}, \code{"p.adjust"}, or \code{"qvalue"}. Default is
#'   \code{"p.adjust"}.
#' @param size Character. Variable for dot size. Must be one of \code{"Count"}
#'   or \code{"GeneRatio"}. Default is \code{"Count"}.
#' @param fontSize Numeric. Base font size for plot text elements.
#'   Default is \code{10}.
#' @param addSeg Logical. Whether to include connecting line segments between
#'   the y-axis and dots (the "stick" of the lollipop). When \code{FALSE},
#'   only dots are displayed. Default is \code{TRUE}.
#' @param lineSize Numeric. Line thickness for the connecting segments when
#'   \code{addSeg = TRUE}. Default is \code{1}.
#' @param lineType Character. Line type for the connecting segments when
#'   \code{addSeg = TRUE}. Options include \code{"solid"}, \code{"dashed"},
#'   \code{"dotted"}, etc. Default is \code{"solid"}.
#'
#' @details
#'   This function creates a lollipop chart where each pathway is represented
#'   by a horizontal line segment (the "stick", controlled by \code{addSeg})
#'   ending in a dot (the "pop"). The visual encoding follows these rules:
#'   \itemize{
#'     \item \strong{X-axis position:} Enrichment metric (NES, Count, or
#'       GeneRatio)
#'     \item \strong{Dot size:} Gene count or gene ratio
#'     \item \strong{Dot/line color:} Enrichment direction (red for enriched,
#'       blue for depleted when x = "NES"; red only for other metrics)
#'     \item \strong{Transparency (alpha):} Statistical significance
#'       (-log10 transformed)
#'     \item \strong{Line segments:} Optional connecting lines from y-axis to
#'       dots, controlled by \code{addSeg}
#'   }
#'
#' @return A named list containing two components:
#'   \describe{
#'     \item{\code{lollipopPlot}}{A \code{ggplot2} object of the lollipop plot
#'       with pathways on the y-axis and enrichment metrics on the x-axis.}
#'     \item{\code{tableLollipopPlot}}{A \code{data.frame} containing the
#'       plotted data with columns:
#'       \itemize{
#'         \item \code{original_name}: Original pathway IDs
#'         \item \code{Description}: Formatted pathway names with line breaks
#'         \item X-axis variable: The column specified by the \code{x} parameter
#'           (\code{NES}, \code{Count}, or \code{GeneRatio})
#'         \item Significance column: The column specified by the \code{alpha}
#'           parameter (\code{pvalue}, \code{p.adjust}, or \code{qvalue})
#'         \item \code{neg_log_alpha}: Transformed significance value
#'           (-log10(alpha))
#'         \item Size variable: The column specified by the \code{size} parameter
#'           (\code{Count} or \code{GeneRatio})
#'       }}
#'   }
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_point scale_color_identity
#'   scale_alpha_continuous scale_size scale_y_discrete scale_x_continuous
#'   guides guide_legend theme_classic theme element_line element_blank
#'   element_text element_rect labs margin
#' @importFrom dplyr filter select arrange
#' @importFrom magrittr %>%
#' @importFrom rlang .data := !!
#' @importFrom methods is
#'
#' @export
#' @examples
#' data("sig2Fun_result")
#' lollipopPlot(seDataFgsea = sig2Fun_result)

lollipopPlot <- function(seDataFgsea, showCategory = 10, breaklineN = 30,
                         x = "NES", alpha = "p.adjust", size = "Count",
                         fontSize = 10, addSeg = TRUE,
                         lineSize = 1, lineType = 'solid') {

    match.arg(x, c("NES", "Count", "GeneRatio"))
    match.arg(alpha, c("pvalue", "p.adjust", "qvalue"))
    match.arg(size, c("Count", "GeneRatio"))

    res.gsea <- .extractDF(seDataFgsea, type = "gseaRaw")
    res <- res.gsea@result

    if (!"Count" %in% colnames(res)) {
        res$Count <- sapply(strsplit(res$core_enrichment, "/"), length)
    }
    if (!"GeneRatio" %in% colnames(res)) {
        res$GeneRatio <- res$Count / res$setSize
    }

    if (is.numeric(showCategory)) {
        if (x == "NES") {
            bottomData <- res %>%
                dplyr::filter(NES <= 0) %>%
                dplyr::arrange(NES) %>%
                dplyr::slice(seq_len(showCategory))

            topData <- res %>%
                dplyr::filter(NES > 0) %>%
                dplyr::arrange(desc(NES)) %>%
                dplyr::slice(seq_len(showCategory))

            res.gsea@result <- rbind(topData, bottomData) %>%
                dplyr::arrange(desc(NES))
        } else {
            res.gsea@result <- res %>%
                dplyr::arrange(desc(!!rlang::sym(x))) %>%
                dplyr::slice(seq_len(showCategory))
        }
    } else {
        res.gsea@result <- res %>%
            dplyr::filter(ID %in% showCategory)
    }

    df <- data.frame(res.gsea)
    data <- if (isS4(res.gsea) && methods::is(res.gsea, "gseaResult")) {
        res.gsea@result
    } else if (is.data.frame(res.gsea)) {
        res.gsea
    } else {
        as.data.frame(res.gsea)
    }

    df <- data.frame(data)

    df$type <- ifelse(df$NES > 0, "activated", "suppressed")

    # Sort data by x-axis variable
    df <- df %>% dplyr::arrange(.data[[x]])
    df$Description <- .labelBreak(df$Description, breaklineN)
    df$Description <- factor(df$Description, levels = df$Description)
    df$type <- factor(df$type, levels = c('suppressed', 'activated'))
    df$neg_log_alpha <- -log10(df[[alpha]])

    # plotting
    if (x == "NES") {
        df$bar_color <- ifelse(df[[x]] >= 0, "#b3200a", "#08306b")
        show_color_legend <- TRUE
    } else {
        df$bar_color <- "#b3200a"
        show_color_legend <- FALSE
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = Description))
    # whether add segment
    if (addSeg == TRUE) {
        p <- p +
            ggplot2::geom_segment(
                ggplot2::aes(x = 0, xend = .data[[x]],
                             y = Description, yend = Description,
                             color = bar_color, alpha = neg_log_alpha),
                linewidth = lineSize, lty = lineType
            )
    }
    p <- p +
        ggplot2::geom_point(
            ggplot2::aes(x = .data[[x]], y = Description,
                         color = bar_color, alpha = neg_log_alpha,
                         size = .data[[size]]),
            stroke = 0
        ) +
        ggplot2::scale_color_identity(
            name = if (show_color_legend) "Enrichment direction" else NULL,
            labels = if (show_color_legend) c("Enriched", "Depleted") else NULL,
            breaks = if (show_color_legend) c("#b3200a", "#08306b") else NULL,
            guide = if (show_color_legend) ggplot2::guide_legend(
                override.aes = list(alpha = 1, size = 5),
                order = 1
            ) else "none",
        ) +
        ggplot2::scale_alpha_continuous(
            name = bquote(-log[10](.(alpha))),
            range = c(0.2, 1),
            labels = function(x) sprintf("%.1f", x),
            guide = ggplot2::guide_legend(
                override.aes = list(color = if (show_color_legend) "black" else "#b3200a", size = 5),
                order = 2
            )
        ) +
        ggplot2::guides(size = ggplot2::guide_legend(
            override.aes = list(shape = 21, fill = "white", color = "black", stroke = 0.6)
        )) +
        ggplot2::scale_size(range = c(3, 7)) +
        ggplot2::labs(x = x, y = NULL) +
        ggplot2::theme_classic(base_size = fontSize) +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_line(color = "grey95"),
            axis.line = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = fontSize, color = "black"),
            axis.text.y = ggplot2::element_text(size = fontSize, color = "black"),
            axis.title.x = ggplot2::element_text(
                colour = "black", size = fontSize + 2, face = "bold",
                margin = margin(t = 5)
            ),
            legend.key = ggplot2::element_rect(fill = "transparent", color = NA),
            legend.background = ggplot2::element_blank(),
            legend.box.background = ggplot2::element_blank()
        ) +
        ggplot2::scale_y_discrete(labels = function(x) .labelBreak(x, breaklineN))

    x_range <- range(df[[x]])
    if (x == "NES") {
        x_buffer <- diff(x_range) * 0.05
        x_min <- min(0, x_range[1] - x_buffer)
        x_max <- max(0, x_range[2] + x_buffer)
        p <- p +
            ggplot2::scale_x_continuous(
                limits = c(x_min, x_max),
                breaks = seq(floor(x_min), ceiling(x_max), by = 1)
            )
    } else {
        x_buffer <- x_range[2] * 0.05
        p <- p +
            ggplot2::scale_x_continuous(
                limits = c(0, x_range[2] + x_buffer),
                expand = c(0, 0)
            )
    }

    tableLollipopPlot <- df |>
        dplyr::select(original_name = ID, Description, !!x, !!alpha, neg_log_alpha, !!size)
    return(list(
        lollipopPlot = p,
        tableLollipopPlot = tableLollipopPlot
    ))
}
