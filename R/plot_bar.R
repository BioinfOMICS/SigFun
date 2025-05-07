#' @title Visualization of Enriched Pathways via Bar Plot
#' @description Generates bar plots displaying significant biological pathways
#' associated with a gene signature, ordered by enrichment significance and
#' colored by association direction. This visualization complements the heatmap
#' output to provide a focused view of top pathway associations.
#'
#' @param SE_data.fgsea A `SummarizedExperiment` object containing:
#' - `fgseaRes` in metadata: GSEA results from `sig2GSEA`
#' - `cor.df` in metadata: Correlation statistics from `sigCor`
#' @param topN Integer. Number of top significant pathways to display per
#' category (default=10)
#' @param significat_type Character. Significance metric for filtering:
#' - "pval" (default): Use nominal p-values
#' - "qval": Use false discovery rate adjusted q-values
#' @param strings Character vector. Pathway categories to visualize
#' (e.g. "HALLMARK", "KEGG"). Default includes major MSigDB collections.
#'
#' @return Returns a list of ggplot objects containing bar plots with:
#' - Pathways ordered by significance (most significant first)
#' - Bar height representing -log10(p/q-value)
#' - Fill color indicating association direction (red = positive, blue = negative)
#' - Horizontal layout optimized for pathway name readability
#'
#' @details Key features:
#' 1. Automatically handles empty results by returning informative empty plots
#' 2. Filters pathways using both statistical significance and string matching
#' 3. Integrates seamlessly with SigFun workflow outputs
#'
#' @examples
#' \dontrun{
#' # After running sig2Fun analysis
#' data("demo_GSE181574")
#' res <- sig2Fun(SE_GSE181574, strings = "HALLMARK")
#'
#' # Generate bar plots for HALLMARK pathways
#' bp <- plot_bar(
#'   SE_data.fgsea = res,
#'   significat_type = "pval",
#'   topN = 10,
#'   strings = "HALLMARK"
#' )
#'
#' # Display plot
#' bp$HALLMARK
#' }
#'
#' @importFrom dplyr filter arrange desc
#' @importFrom S4Vectors metadata
#' @importFrom ggplot2 ggplot theme_void annotate lims theme unit
#' @export

plot_bar <- function(SE_data.fgsea, significat_type="pval",
    topN=10, strings=c("GOBP", "GOCC", "GOMF", "KEGG", "REACTOME", "WP")) {

    RES_GSEA <- S4Vectors::metadata(SE_data.fgsea)$fgsea %>%
                    dplyr::filter(pval < 0.05)

    pattern <- paste0("[", paste(paste(strings, " "), collapse="|"), "]")

    RES_NES_total <- RES_GSEA %>% dplyr::slice(grep(pattern, pathway)) %>%
        dplyr::mutate(`-log10(pvalue)`=-log10(pval))

    RES_NES_total$leadingEdge <- gsub("[|]", ";", RES_NES_total$leadingEdge)

    #data.table::fwrite(RES_NES_total, col.names=TRUE, row.names=FALSE,
    #file=file.path(output_path, "GSEA_stat.txt"), sep="\t", quote=FALSE)
    barplots <- NULL
    for (j in seq_along(strings)) {

        #plot.name <- paste0("GSEA_top10_", gsub(" ", "", strings[j]), ".png")

        RES_NES_strings <- RES_NES_total %>%
            dplyr::filter(grepl(paste0(strings[j], " "), pathway),
                        .data[[significat_type]] < 0.05) %>%
            dplyr::arrange(desc(NES))

        if (nrow(RES_NES_strings) > 0) {

        barplots[[strings[j]]] <- .barplot(plot.name=plot.name,
                    type.sig=significat_type, topN=10,
                    RES_NES_strings=RES_NES_strings)

        } else {

            p <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate(
    geom="text",x=1, y=1, label=paste("No significant (",
                significat_type, "<0.05) item.")) +
                ggplot2::lims(x=c(0, 2), y=c(0, 2)) +
                ggplot2::theme(plot.margin=ggplot2::unit(c(0, 0, 0, 0), "cm"))
            barplots[[strings[j]]] <- p
            #ggplot2::ggsave(filename=file.path(output_path, plot.name),
            #    plot=p,width=3,height=1.4,units="in",dpi=600)
        }
    }
    return(barplots)
}


