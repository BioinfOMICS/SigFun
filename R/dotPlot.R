#' @title dotPlot
#' @description The \code{dotPlot} function generates a customizable dot plot
#'   visualization for Gene Set Enrichment Analysis (GSEA) results stored in a
#'   \code{SummarizedExperiment} object. The plot displays enriched pathways as
#'   dots positioned by enrichment metrics (e.g., \code{GeneRatio}, \code{Count},
#'   or \code{NES}) and encoded by color, size, and transparency to represent
#'   statistical significance and biological importance.
#'
#' @param seDataFgsea A \code{SummarizedExperiment} object containing GSEA
#'   results. The metadata must include a \code{gseaResult} object with pathway
#'   information, gene sets, enrichment statistics, and ranking metrics.
#' @param showCategory Either a numeric value specifying the number of top
#'   enriched pathways to display, or a character vector specifying particular
#'   pathway IDs or names. Default is \code{10}.
#' @param breaklineN Integer. Maximum number of characters before inserting
#'   line breaks into pathway names for better readability. Default is \code{30}.
#' @param x Character. Variable to display on the x-axis. Must be one of
#'   \code{"GeneRatio"}, \code{"Count"}, or \code{"NES"}. Default is
#'   \code{"GeneRatio"}.
#' @param size Character. Variable used to scale dot size. Must be one of
#'   \code{"Count"} or \code{"GeneRatio"}. Default is \code{"Count"}.
#' @param alpha Character. Statistical metric used for transparency encoding.
#'   Must be one of \code{"pvalue"}, \code{"p.adjust"}, or \code{"qvalue"}.
#'   Default is \code{"p.adjust"}.
#' @param orderBy Character. Column name used to order pathways before plotting.
#'   Default is \code{"GeneRatio"}.
#' @param decreasing Logical. Whether to sort pathways in decreasing order.
#'   Default is \code{TRUE}.
#' @param fontSize Numeric. Base font size for all plot text elements.
#'   Default is \code{10}.
#' @param title Character or \code{NULL}. Custom plot title. If \code{NULL},
#'   no title is displayed. Default is \code{NULL}.
#'
#' @details
#'   This function extends the visualization of enrichment results beyond the
#'   default \pkg{enrichplot::dotplot()} by supporting flexible mapping of dot
#'   color, size, and alpha transparency to different variables.
#'
#'   Dots represent enriched gene sets, where:
#'   \itemize{
#'     \item \code{x}-axis: Enrichment metric (e.g., GeneRatio, Count, NES)
#'     \item \code{Dot size}: Pathway size or gene ratio
#'     \item \code{Dot color}: Enrichment direction (for NES) or fixed color
#'     \item \code{Transparency (alpha)}: Statistical significance
#'   }
#'
#'   Pathways can be selected either by the number of top categories
#'   (\code{showCategory = 10}) or explicitly by ID/name.
#'   When \code{x = "NES"}, the plot uses a diverging red–blue palette to
#'   indicate up- and down-regulated enrichment directions.
#'
#' @return A named list with two components:
#'   \describe{
#'     \item{\code{dotPlot}}{A \code{ggplot2} object visualizing enriched pathways
#'       with customizable color, size, and transparency encodings.}
#'     \item{\code{tableDotPlot}}{A \code{data.frame} containing the plotted data,
#'       including columns:
#'       \itemize{
#'         \item \code{original_name}: Original pathway IDs
#'         \item \code{Description}: Formatted pathway names (with line breaks)
#'         \item \code{NES}: Normalized enrichment score
#'         \item \code{Count}: Number of core genes per pathway
#'         \item \code{GeneRatio}: Ratio of core genes to set size
#'         \item Statistical significance column: The column specified by the
#'           \code{alpha} parameter (\code{pvalue}, \code{p.adjust}, or \code{qvalue})
#'       }}
#'   }
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_identity
#'   scale_alpha_continuous scale_size scale_y_discrete theme_minimal theme
#'   labs guides guide_legend element_blank element_text element_line
#'   element_rect margin
#' @importFrom dplyr mutate select all_of
#' @importFrom scales rescale
#' @importFrom rlang .data
#'
#' @export
#' @examples
#' data("sig2Fun_result")
#' dotPlot(seDataFgsea = sig2Fun_result)
#'
#' @note
#' Requires \pkg{ggplot2} (>= 3.5.0) for legend customization and margin support.

# Note: Requires ggplot2 (>= 3.5.0)
dotPlot <- function(
        seDataFgsea, showCategory = 10, breaklineN = 30,
        x = "GeneRatio", size = "Count", alpha = "p.adjust",
        orderBy = "GeneRatio", decreasing = TRUE,
        fontSize = 10, title = NULL
) {

    match.arg(alpha, c("pvalue", "p.adjust", "qvalue"))
    match.arg(x, c("GeneRatio", "Count", "NES"))
    match.arg(size, c("Count", "GeneRatio"))

    gseaRaw <- .extractDF(seDataFgsea, type = "gseaRaw")
    df <- gseaRaw@result |>
        dplyr::mutate(Description = .labelBreak(Description, breaklineN))

    if (!("Count" %in% colnames(df)))
        df$Count <- sapply(strsplit(df$core_enrichment, "/"), length)
    if (!("GeneRatio" %in% colnames(df)))
        df$GeneRatio <- df$Count / df$setSize

    if (is.numeric(showCategory)) {
        df <- df[seq_len(min(showCategory, nrow(df))), ]
    } else {
        df <- df[df$ID %in% showCategory, ]
    }

    if (nrow(df) == 0) {
        cli::cli_abort("No pathways to display. Please check {.arg showCategory}.")
    }

    if (orderBy %in% colnames(df)) {
        df <- df[order(df[[orderBy]], decreasing = decreasing), ]
    }

    df$Description <- factor(df$Description, levels = rev(unique(df$Description)))
    df$neg_log_alpha <- -log10(df[[alpha]])
    max_neg_log <- max(df$neg_log_alpha, na.rm = TRUE)
    if (!is.finite(max_neg_log) || max_neg_log <= 0) {
        max_neg_log <- 1
    }
    df$alpha_value <- sqrt(df$neg_log_alpha / max_neg_log)
    df$alpha_value <- scales::rescale(df$alpha_value, to = c(0.2, 1))

    if (x == "NES") {
        df$dot_color <- ifelse(df$NES >= 0, "#b3200a", "#08306b")
    } else {
        df$dot_color <- "#b3200a"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[["Description"]])) +
        ggplot2::geom_point(
            ggplot2::aes(size = .data[[size]], color = dot_color, alpha = neg_log_alpha),
            stroke = 0
        ) +
        ggplot2::scale_color_identity(
            name = if (x == "NES") "Enrichment direction" else NULL,
            labels = if (x == "NES") c("Enriched", "Depleted") else NULL,
            breaks = if (x == "NES") c("#b3200a", "#08306b") else NULL,
            guide = if (x == "NES") ggplot2::guide_legend(
                override.aes = list(alpha = 1, size = 6), order = 1)
            else "none"
        ) +
        ggplot2::scale_alpha_continuous(
            name = bquote(-log[10](.(alpha))),
            range = c(0.2, 1),
            labels = function(x) sprintf("%.1f", x),
            guide = ggplot2::guide_legend(
                override.aes = list(color = if (x == "NES") "black" else "#b3200a", size = 6),
                order = 2
            )
        ) +
        ggplot2::scale_size(range = c(3, 8)) +
        ggplot2::scale_y_discrete(labels = function(x) .labelBreak(x, breaklineN)) +
        ggplot2::theme_minimal(base_size = fontSize) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "grey95", linewidth = 0.3),
            plot.title = ggplot2::element_text(size = fontSize + 4, face = "bold",
                                               hjust = 0.5, margin = ggplot2::margin(b = 10)),
            axis.text.y = ggplot2::element_text(size = fontSize, color = "black"),
            axis.title.x = ggplot2::element_text(size = fontSize + 2, face = "bold",
                                                 margin = ggplot2::margin(t = 8)),
            legend.key = ggplot2::element_rect(fill = "transparent", color = NA),
            legend.background = ggplot2::element_blank(),
            legend.box.background = ggplot2::element_blank()
        ) +
        ggplot2::labs(x = x, y = NULL, title = title) +
        ggplot2::guides(
            size = ggplot2::guide_legend(
                override.aes = list(shape = 21, fill = "white", color = "black", stroke = 0.6))
        )

    tableDotPlot <- df |>
        dplyr::select(original_name = ID, Description, NES, Count, GeneRatio, dplyr::all_of(alpha))

    return(list(dotPlot = p, tableDotPlot = tableDotPlot))
}
