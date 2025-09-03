#' @title barplot
#' @description The \code{barplot} function generates a horizontal bar plot
#' for visualizing gene set enrichment analysis (GSEA) results. It displays the
#' gene count for each pathway on the x-axis and pathway names on the y-axis,
#' with bars colored by statistical significance measures or enrichment scores.
#' @param SE_data.fgsea A \code{SummarizedExperiment} object containing GSEA
#'   results. The object must contain:
#'   \itemize{
#'     \item \code{metadata(SE_data.fgsea)$gseaResult}: GSEA result object with
#'           result slot containing pathway information
#'   }
#' @param topN Integer. Number of top significant pathways to display per
#' category (default=10).
#' @param breaklineN An integer specifying the number of characters after which
#' to break pathway names into multiple lines for better readability on the
#' y-axis. Default is \code{30}.
#' @param alpha A character string specifying which statistical measure to use
#' for bar coloring transparency. Must be one of
#' \code{c('pvalue', 'p.adjust', 'NES')}. Default uses argument matching.
#'   \itemize{
#'     \item \code{pvalue}: Raw p-values
#'     \item \code{p.adjust}: Adjusted p-values (multiple testing correction)
#'     \item \code{NES}: Normalized Enrichment Score
#'   }
#' @return A named list containing two components:
#' \itemize{
#'   \item \code{barPlot}: A \code{ggplot2} object representing a horizontal
#'         bar plot with:
#'         \itemize{
#'           \item X-axis: Gene count per pathway
#'           \item Y-axis: Formatted pathway descriptions (ordered by count)
#'           \item Bar color: Selected statistical measure (p-value, adjusted
#'                 p-value, or NES)
#'           \item Color legend: Continuous scale for the selected measure
#'         }
#'   \item \code{Table_barPlot}: A \code{data.frame} containing the underlying
#'         data used for plotting with columns:
#'         \itemize{
#'           \item \code{original_name}: Original pathway identifiers
#'           \item \code{label_name}: Formatted pathway names with line breaks
#'           \item \code{Count}: Gene count for each pathway
#'           \item Column named by \code{color} parameter: Statistical measure
#'                 values used for coloring
#'         }
#' }
#' @export
#' @examples
#' data("demo_GSE181574")
#' barPlot(SE_data.fgsea=GSE181574.sigfun)
barPlot <- function(
    SE_data.fgsea, topN=10, breaklineN=30, alpha='pvalue'){
    match.arg(alpha, c('pvalue', 'p.adjust', 'qvalue'))
    gseaRaw <- .extractDF(SE_data.fgsea, type='gseaRaw')
    bottom.NES <- gseaRaw@result %>%
      dplyr::filter(NES <= 0) %>%
      dplyr::arrange(NES) %>% dplyr::slice(seq_len(topN))
    top.NES <- gseaRaw@result %>%
      dplyr::filter(NES > 0) %>%
      dplyr::arrange(desc(NES)) %>% dplyr::slice(seq_len(topN))
    all.NES <- rbind(top.NES, bottom.NES) %>%
      dplyr::arrange(desc(NES)) %>%
      dplyr::select(ID, NES, dplyr::all_of(alpha))
    all.NES$alpha <- all.NES[,3]
    all.NES <- all.NES%>%
      dplyr::mutate(name=.labelBreak(ID, breaklineN)) %>%
      dplyr::mutate(name=forcats::fct_reorder(name, NES)) %>%
      dplyr::mutate(
        neg_log_alpha=-log10(alpha),
        neg_log_alpha=ifelse(is.infinite(neg_log_alpha),
          max(neg_log_alpha[is.finite(neg_log_alpha)], na.rm=TRUE),
          neg_log_alpha),
        bar_color=ifelse(alpha >= 0.05, 'grey60',
          ifelse(NES >= 0, '#FF5151', '#4169E1')))
    max_neg_log <- max(all.NES$neg_log_alpha, na.rm=TRUE)
    all.NES$alpha_value <- pmax(0, all.NES$neg_log_alpha / max_neg_log)
    barPlot <- ggplot2::ggplot(all.NES,
      ggplot2::aes(y=name, x=NES, fill=bar_color, alpha=alpha_value)) +
      geom_col(width=.8) +
      ggplot2::scale_fill_identity(
        name="Significance", labels=c("Enriched", "NS", "Depleted"),
        breaks=c('#FF5151', "grey60", '#4169E1'),
        guide=ggplot2::guide_legend(override.aes=list(alpha=1))) +
      ggplot2::scale_alpha_identity(
        name=paste0("-log10(", alpha, ")"),
        breaks=seq(0, 1, 0.25),
        labels=sprintf("%.1f", seq(0, 1, 0.25) * max_neg_log),
        guide=ggplot2::guide_legend(override.aes=list(fill="black"))) +
      ggplot2::scale_x_continuous(
        expand=ggplot2::expansion(mult=c(0.02, 0.02))) +
      ggplot2::labs(x=NULL, Y=NULL, title=NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid=ggplot2::element_blank(),
        axis.title.x=ggplot2::element_text(
          size=12, margin=ggplot2::margin(t=10)),
        axis.title.y=ggplot2::element_blank(),
        axis.text.y=ggplot2::element_text(
          color="black", size=8, hjust=1, lineheight=.8),
        axis.text.x=ggplot2::element_text(size=10),
        axis.line.x=ggplot2::element_line(color="black", linewidth=0.5),
        axis.line.y=ggplot2::element_line(color="black", linewidth=0.5),
        axis.ticks.x=ggplot2::element_line(color="black", linewidth=0.3),
        axis.ticks.y=ggplot2::element_line(color="grey40", linewidth=0.5),
        axis.ticks.length.y=grid::unit(0.2, "cm"),
        axis.ticks.length.x=grid::unit(0.15, "cm"),
        legend.position="right",
        legend.title=ggplot2::element_text(size=10),
        legend.text=ggplot2::element_text(size=8),
        legend.key.height=grid::unit(0.8, "cm"),
        legend.key.width=grid::unit(0.8, "cm"),
        legend.spacing.y=grid::unit(0.5, "cm"),
        plot.margin=ggplot2::margin(10, 10, 10, 10)
      )
    Table_barPlot <- all.NES %>%
      dplyr::select(original_name=ID, label_name=name, dplyr::everything())
    return(list(barPlot=barPlot, Table_barPlot=Table_barPlot))
}
