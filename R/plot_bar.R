#' @title plot_bar
#' @description This function is explore the functions related to signature by a
#' surrogate of whole transcriptome data.
#' @param SE_data.fgsea Summarized Experiment Object. In the metadata of the
#' input, it must be a Summary Experiment Object and should contain the
#' following two objects "cor.df" and "fgseaRes".
#' "cor.df" is the result of SigCor. Add the cor.df to the metadata by running
#' SigCor. "fgseaRes" is the result of Sig2GSEA.
#' Add the fgseaRes to the metadata by running Sig2GSEA
#' @param output_path Character. Location of the output path.
#' @param topN Integer. Number of significant functions in hte output bar plots.
#' Default is 10.
#' @param significat_type Character. Filter or statistical significance.
#' "pval" for p-value and "qval" for Q-value
#' @param strings vector. Pathways header want to be plotted.
#' @return Return bar plots to the outpath specified.
#' @export
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @examples
#' data("demo")
#' output_path <- tempdir()
#' plot_bar(
#' SE_data.fgsea = SE_data,
#' output_path=output_path,
#' topN=10,
#' significat_type = "pval",
#' strings = c("KEGG"))
plot_bar <- function(SE_data.fgsea, output_path, significat_type = "pval",
    topN = 10, strings = c("GOBP", "GOCC", "GOMF", "KEGG", "REACTOME", "WP")) {

    RES_GSEA <- SE_data.fgsea@metadata$fgseaRes %>%
                    dplyr::filter(pval < 0.05)

    pattern <- paste0("[", paste(paste(strings, " "), collapse = "|"), "]")

    RES_NES_total <- RES_GSEA %>% dplyr::slice(grep(pattern, pathway)) %>%
        dplyr::mutate(`-log10(pvalue)` = -log10(pval))

    RES_NES_total$leadingEdge <- gsub("[|]", ";", RES_NES_total$leadingEdge)

    data.table::fwrite(RES_NES_total, col.names = TRUE, row.names = FALSE,
    file = file.path(output_path, "GSEA_stat.txt"), sep = "\t", quote = FALSE)
    barplots <- NULL
    for (j in seq_along(strings)) {

        plot.name <- paste0("GSEA_top10_", gsub(" ", "", strings[j]), ".png")

        RES_NES_strings <- RES_NES_total %>%
            dplyr::filter(grepl(paste0(strings[j], " "), pathway),
                        .data[[significat_type]] < 0.05) %>%
            dplyr::arrange(desc(NES))

        if (nrow(RES_NES_strings) > 0) {

        barplots[[strings[j]]] <- .barplot(plot.name=plot.name,
                    type.sig=significat_type, topN=10,
                    RES_NES_strings=RES_NES_strings, output_path=output_path)

        } else {

            p <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate(
    geom = "text",x = 1, y = 1, label = paste("No significant (",
                significat_type, "<0.05) item.")) +
                ggplot2::lims(x = c(0, 2), y = c(0, 2)) +
                ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"))

            ggplot2::ggsave(filename = file.path(output_path, plot.name),
                plot = p,width = 3,height = 1.4,units = "in",dpi = 600)
        }
    }
    return(barplots)
}


