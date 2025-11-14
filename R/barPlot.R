#' @title barPlot
#' @description The \code{barPlot} function generates a horizontal bar plot
#'   for visualizing gene set enrichment analysis (GSEA) results. It supports
#'   multiple x-axis variables (e.g., NES, Count, or GeneRatio) and allows
#'   bar coloring and transparency to be mapped by statistical significance
#'   measures (p-value, adjusted p-value, or q-value). This function provides
#'   a flexible and publication-ready visualization of pathway-level enrichment.
#'
#' @param seDataFgsea A \code{SummarizedExperiment} object containing GSEA
#'   results. The object must contain:
#'   \itemize{
#'     \item \code{metadata(seDataFgsea)$gseaResult}: GSEA result object whose
#'       result slot contains pathway-level enrichment information.
#'   }
#' @param topN Numeric. Number of top-ranked pathways to display.
#'   Default is \code{10}.
#' @param breaklineN Integer. Maximum number of characters before wrapping
#'   pathway names for improved readability. Default is \code{30}.
#' @param alpha Character. Variable for transparency encoding.
#'   Must be one of \code{'pvalue'}, \code{'p.adjust'}, or \code{'qvalue'}.
#'   \itemize{
#'     \item \code{pvalue}: Raw p-values.
#'     \item \code{p.adjust}: Adjusted p-values (multiple testing correction).
#'     \item \code{qvalue}: q-values (alternative adjusted significance measure).
#'   }
#' @param x Character. Variable for the x-axis.
#'   Must be one of \code{'NES'}, \code{'Count'}, or \code{'GeneRatio'}.
#'   \itemize{
#'     \item \code{NES}: Normalized Enrichment Score (default).
#'     \item \code{Count}: Number of core enrichment genes in each pathway.
#'     \item \code{GeneRatio}: Ratio of enriched genes to total genes in the set.
#'   }
#' @param fontSize Numeric. Base font size for plot text elements.
#'   Default is \code{8}.
#'
#' @importFrom magrittr %>%
#' @importFrom forcats fct_reorder
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_identity scale_alpha_continuous
#'   scale_x_continuous labs theme_classic theme element_text element_line element_blank
#'   guide_legend expansion margin
#' @importFrom dplyr filter arrange mutate select slice all_of
#'
#' @return A named list containing two components:
#'   \describe{
#'     \item{\code{barPlot}}{A \code{ggplot2} object representing the horizontal
#'       bar plot, where:
#'       \itemize{
#'         \item X-axis: Selected variable (\code{NES}, \code{Count}, or
#'           \code{GeneRatio})
#'         \item Y-axis: Formatted pathway names (ordered by the selected variable)
#'         \item Bar color: Encodes enrichment direction or significance
#'         \item Transparency (alpha): Proportional to \code{-log10(alpha)} of the
#'           chosen statistic
#'       }}
#'     \item{\code{tableBarPlot}}{A \code{data.frame} containing the underlying
#'       plotting data with columns:
#'       \itemize{
#'         \item \code{original_name}: Original pathway identifiers.
#'         \item \code{label_name}: Formatted pathway names with line breaks.
#'         \item \code{NES / Count / GeneRatio}: The variable used on the x-axis,
#'           depending on the selected \code{x} parameter.
#'         \item \code{alpha}: Statistical values used for transparency
#'           (e.g., p-value, adjusted p-value, or q-value).
#'         \item \code{neg_log_alpha}: The transformed significance level
#'           calculated as \code{-log10(alpha)}.
#'         \item \code{alpha_value}: Rescaled transparency value ranging from 0–1,
#'           used for plotting alpha intensity.
#'       }}
#'   }
#'
#' @export
#' @examples
#' data("sig2Fun_result")
#' barPlot(seDataFgsea = sig2Fun_result)
#'
#' @note
#' Requires \pkg{ggplot2} (>= 3.5.0) and \pkg{dplyr} (>= 1.1.0) for full compatibility.

# Note: Requires ggplot2 (>= 3.5.0) and dplyr (>= 1.1.0)
barPlot <- function(seDataFgsea, topN = 10, breaklineN = 30, alpha = "pvalue",
                    x = "NES", fontSize = 8) {

    match.arg(alpha, c("pvalue", "p.adjust", "qvalue"))
    match.arg(x, c("NES", "Count", "GeneRatio"))

    gseaRaw <- .extractDF(seDataFgsea, type = "gseaRaw")
    res <- gseaRaw@result
    if (!("Count" %in% colnames(res))) {
        res$Count <- sapply(strsplit(res$core_enrichment, "/"), length)
    }
    if (!("GeneRatio" %in% colnames(res))) {
        res$GeneRatio <- res$Count / res$setSize
    }

    if (x == "NES") {
        bottomData <- res %>%
            dplyr::filter(NES <= 0) %>%
            dplyr::arrange(NES) %>% dplyr::slice(seq_len(topN))
        topData <- res %>%
            dplyr::filter(NES > 0) %>%
            dplyr::arrange(desc(NES)) %>% dplyr::slice(seq_len(topN))
        allData <- rbind(topData, bottomData) %>%
            dplyr::arrange(desc(NES))
    } else {
        allData <- res %>%
            dplyr::arrange(desc(!!rlang::sym(x))) %>%
            dplyr::slice(seq_len(topN))
    }

    allData <- allData %>% dplyr::select(ID, dplyr::all_of(x), dplyr::all_of(alpha))
    allData$alpha <- allData[[alpha]]
    allData[[alpha]][is.na(allData[[alpha]])] <- 1
    allData <- allData %>%
        dplyr::mutate(
            name = .labelBreak(ID, breaklineN),
            name = forcats::fct_reorder(name, !!rlang::sym(x)),
            neg_log_alpha = -log10(alpha),
            neg_log_alpha = ifelse(is.infinite(neg_log_alpha),
                                   max(neg_log_alpha[is.finite(neg_log_alpha)], na.rm = TRUE),
                                   neg_log_alpha)
        )

    if (x == "NES") {
        allData <- allData %>%
            dplyr::mutate(
                bar_color = ifelse(alpha >= 0.05, "grey60",
                                   ifelse(NES >= 0, "#D25C43", "#5979A3"))
            )
    } else {
        allData <- allData %>%
            dplyr::mutate(
                bar_color = ifelse(alpha >= 0.05, "grey60", "#D25C43")
            )
    }

    max_neg_log <- max(allData$neg_log_alpha, na.rm = TRUE)
    if (is.na(max_neg_log) || is.infinite(max_neg_log) || max_neg_log <= 0) {
        message("No significant pathway detected — using default alpha scale.")
        max_neg_log <- 1
    }
    allData$alpha_value <- pmax(0, allData$neg_log_alpha / max_neg_log)

    bar_width <- if (x == "NES") 0.7 else 0.6
    barPlot <- ggplot2::ggplot(
        allData,
        ggplot2::aes(
            y = name, x = !!rlang::sym(x),
            fill = bar_color, alpha = neg_log_alpha
        )
    ) +
        ggplot2::geom_col(width = bar_width) +
        ggplot2::scale_fill_identity(
            name = "Significance",
            labels = if (x == "NES") c("Enriched", "NS", "Depleted") else c("Significant", "NS"),
            breaks = if (x == "NES") c("#D25C43", "grey60", "#5979A3") else c("#D25C43", "grey60"),
            guide = ggplot2::guide_legend(override.aes = list(alpha = 1))
        ) +
        ggplot2::scale_alpha_continuous(
            name = bquote(-log[10](.(alpha))),
            range = c(0.2, 1),
            guide = ggplot2::guide_legend(override.aes = list(fill = "black"))
        ) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
        ggplot2::labs(x = x, y = NULL, title = NULL) +
        ggplot2::theme_classic(base_size = fontSize) +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_text(size = fontSize + 2, face = "bold",
                                                 margin = ggplot2::margin(t = 5)),
            axis.text.y = ggplot2::element_text(color = "black", size = fontSize, hjust = 1),
            axis.text.x = ggplot2::element_text(color = "black", size = fontSize, hjust = 0.5),
            axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.6),
            axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.6),
            axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
            legend.position = "right",
            legend.title = ggplot2::element_text(size = fontSize),
            legend.text = ggplot2::element_text(size = fontSize),
            legend.key.height = grid::unit(0.5, "cm"),
            legend.key.width = grid::unit(0.5, "cm"),
            plot.margin = ggplot2::margin(5, 5, 5, 5)
        )

    tableBarPlot <- allData %>%
        dplyr::select(original_name = ID, label_name = name, !!x, !!alpha, neg_log_alpha, alpha_value)
    return(list(barPlot = barPlot, tableBarPlot = tableBarPlot))
}
