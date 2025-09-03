#' @title heatPlot
#' @description The \code{heatplot} function generates a dot plot heatmap for
#' visualizing gene set enrichment analysis (GSEA) results. It displays genes
#' on the x-axis and pathways on the y-axis, with dots colored by correlation
#' coefficients and sized by statistical significance measures.
#' @param SE_data.fgsea A \code{SummarizedExperiment} object containing GSEA
#'   results and correlation data. The object must contain:
#'   \itemize{
#'     \item \code{metadata(SE_data.fgsea)$gseaResult}: GSEA result object
#'     \item \code{metadata(SE_data.fgsea)$cor.df}: Correlation data frame with
#'           gene, correlation, and statistical columns
#'   }
#' @param breaklineN An integer specifying the number of characters after which
#'   to break pathway names into multiple lines for better readability.
#'   Default is \code{30}.
#' @param showCategory Either a numeric value indicating the number of top
#'   pathways to display, or a character vector specifying particular pathway
#'   names to show. Default is \code{5}.
#' @param distfun Character. The distance measure to be used when choosing
#' \code{"hclustering"} as clustering method. Allow method include "euclidean",
#' "manhattan", "maximum", "canberra", "binary", and "minkowski".
#' @param hclustfun Character. The agglomeration method to be used when choosing
#' \code{"hclustering"} as clustering method. This should be (an unambiguous
#' abbreviation of) one of "ward.D", "ward.D2", "single", "complete",
#' "average" (=UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC), or "centroid"
#' (= UPGMC).
#' @param dotColor A character string specifying which column from the
#' correlation data frame to use for dot color.
#' @param dotSize A character string specifying which column from the
#' correlation data frame to use for dot sizing. Must be one of
#' \code{c('pval', 'p.adjust')}. Default is \code{pval}.
#' @return A list containing two elements:
#'   \describe{
#'     \item{\code{heatPlot}}{A ggplot2 object representing a dot plot heatmap
#'           with genes clustered on the x-axis and pathways with formatted
#'           names on the y-axis. Dots are colored by correlation coefficients
#'           and sized by statistical significance measures}
#'     \item{\code{table_heatPlot}}{A data.frame containing the plot data
#'           with columns: \code{original_name} (original pathway IDs),
#'           \code{label_name} (formatted pathway names with line breaks),
#'           \code{Gene} (gene symbols), \code{Coefficient} (correlation
#'           coefficients), and a column named after the \code{dotSize}
#'           parameter (statistical significance values used for dot sizing)}
#'   }
#' @export
#' @examples
#' data("demo_GSE181574")
#' heatPlot(SE_data.fgsea=GSE181574.sigfun, showCategory=3)
heatPlot <- function(SE_data.fgsea, breaklineN=30, showCategory=5,
    dotSize='pval', dotColor='cor', distfun='euclidean', hclustfun='ward.D'){
    match.arg(dotSize, c('pval', 'p.adjust'))
    gseaReadable <- .extractDF(SE_data.fgsea, type="gseaReadable")
    cor.df <- .extractDF(SE_data.fgsea, type="corCoef")
    if(isFALSE(dotColor %in% colnames(cor.df))){
        cli::cli_abort(
        'The input {.arg dotColor} must be a {.cls character} and must be
        included in the column names of {.field SE_data.fgsea@metadata$cor.df}.')
    }
    if(isFALSE(dotSize %in% colnames(cor.df))){
      cli::cli_abort(
        'The input {.arg dotSize} must be a {.cls character} and must be
        included in the column names of {.field SE_data.fgsea@metadata$cor.df}.')
    }
    data <- .extractGeneSets(gseaReadable, showCategory)
    fontSize <- ifelse(length(unique(data$categoryID)) < 10, 12, 10)
    cor.df <- cor.df |>
      dplyr::select(gene, dplyr::all_of(dotColor), dplyr::all_of(dotSize)) |>
      dplyr::rename(Gene=gene,
        Coef=dplyr::all_of(dotColor), size=dplyr::all_of(dotSize))
    if (gseaReadable@readable && gseaReadable@keytype != "SYMBOL") {
        cor.df <- data.frame(Gene=names(gseaReadable@gene2Symbol),
            genesymbol=gseaReadable@gene2Symbol) |>
            merge(cor.df) |> dplyr::select(-Gene) |>
            dplyr::rename(Gene=genesymbol)
    }
    data <- dplyr::left_join(data, cor.df) |>
        dplyr::mutate(name=.labelBreak(categoryID, breaklineN))
    heatPlot <- data |>
        dplyr::mutate(Gene=factor(Gene,
            levels=.geneOrder(data, distfun, hclustfun))) |>
        ggplot2::ggplot(ggplot2::aes(x=Gene, y=name)) +
        ggplot2::geom_point(ggplot2::aes(size=size, fill=Coef),
            color='black', shape=21) +
        .scale_fill(data$Coef) + ggplot2::xlab(NULL) + ggplot2::ylab(NULL) +
        ggplot2::guides(fill=ggplot2::guide_colorbar(title=dotColor)) +
        ggplot2::theme_minimal() +
        ggplot2::scale_size(name=dotSize, range=c(4, 1.2)) +
        ggplot2::theme(
            panel.grid.major=ggplot2::element_blank(),
            axis.text.x=ggplot2::element_text(
                angle=60, hjust=1, color='black', size=fontSize),
            axis.text.y=ggplot2::element_text(color='black', size=12))
    table_heatPlot <- data |>
        dplyr::select(original_name=categoryID, label_name=name,
            Gene, Coefficient=Coef, !!dotSize := size)
    return(list(heatPlot=heatPlot, table_heatPlot=table_heatPlot))
}
