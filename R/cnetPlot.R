#' @title cnetPlot
#' @description Generates a concept network visualization for gene set
#' enrichment analysis(GSEA) results. The network displays connections between
#' enriched pathways and their associated genes, with optional color coding
#' based on gene expression correlation values. This visualization helps
#' identify hub genes that participate in multiple biological pathways and
#' reveals the interconnected structure of enriched gene sets.
#' @param SE_data.fgsea A \code{SummarizedExperiment} object containing GSEA
#'   results in its metadata. Required metadata components:
#'   \describe{
#'     \item{\code{gseaResult}}{GSEA result object with pathway information,
#'       gene sets, and enrichment statistics}
#'     \item{\code{cor.df}}{Data frame with gene correlation values (when
#'       \code{color=TRUE}), containing \code{gene} and \code{cor} columns}
#'   }
#' @param showCategory Either a numeric value indicating the number of top
#' pathways to display, or a character vector specifying particular pathway
#' names to show. Default is \code{5}.
#' @param node_label Character string specifying node label type. One of 'all',
#' 'none', 'category', 'item', 'exclusive' or 'share'.
#' Default: \code{"category"}
#' @param fontSize Numeric value controlling category label font size in the
#' network. Default: \code{1.2}
#' @param color Logical indicating whether to color-code genes by their
#' correlation values. Requires correlation data in metadata when \code{TRUE}.
#' Default: \code{TRUE}
#' @param dotColor A character string specifying which column from the
#' correlation data frame to use for dot color.
#' @param circular Logical indicating whether to use circular layout
#' (\code{TRUE}) or automatic nice layout (\code{FALSE}). Default: \code{FALSE}
#' @return A list containing three elements:
#'   \describe{
#'     \item{\code{cnetPlot}}{A ggplot2 object representing the concept network
#'           visualization with pathways and genes connected by edges,
#'           optionally colored by correlation values}
#'     \item{\code{TableEdge_cnetPlot}}{A data.frame containing edge information
#'           with columns: \code{from} (pathway IDs) and \code{to} (individual
#'           genes from core enrichment), representing the connections between
#'           pathways and their associated genes}
#'     \item{\code{TableNode_cnetPlot}}{A data.frame containing node information
#'           with columns: \code{id} (node names/identifiers), \code{size}
#'           (node size values), and optionally \code{foldChange} (correlation
#'           values for color coding, included only when \code{color=TRUE})}
#'   }
#' @param breaklineN An integer specifying the number of characters after which
#' to break pathway names into multiple lines for better readability.
#' @export
#' @examples
#' data("demo_GSE181574")
#' cnetPlot(SE_data.fgsea=GSE181574.sigfun)
cnetPlot <- function(SE_data.fgsea, showCategory=5, node_label='category',
    fontSize=1.2, color=TRUE, dotColor='cor', circular=FALSE, breaklineN=30){
    res.gsea <- .extractDF(SE_data.fgsea, type='gseaReadable')
    if(color){
        cor.df <- .extractDF(SE_data.fgsea, type='corCoef')
        if (res.gsea@readable && res.gsea@keytype != "SYMBOL") {
          cor.df <- data.frame(gene=names(res.gsea@gene2Symbol),
            genesymbol=res.gsea@gene2Symbol) |>
            merge(cor.df) |> dplyr::select(-gene) |>
            dplyr::rename(gene=genesymbol)
        }
        if(isFALSE(dotColor %in% colnames(cor.df))){
          cli::cli_abort(
            'The input {.arg dotColor} must be a {.cls character} and must be
        included in the column names of {.field SE_data.fgsea@metadata$cor.df}.')
        }else{
          geneList <- cor.df |> dplyr::pull(dotColor) |> setNames(cor.df$gene) |>
            sort(decreasing=TRUE)
        }
    }else{
        geneList <- NULL
    }
    #if (is.numeric(showCategory)) {
    #    res.gsea@result <- res.gsea@result[seq_len(showCategory), ]
    #}else{
    #    res.gsea@result <- res.gsea@result |>
    #        dplyr::filter(ID %in% showCategory)
    #}
    geneSets <- DOSE::geneInCategory(res.gsea)[res.gsea@result$ID]
    if (is.numeric(showCategory)){
      geneSets <- geneSets[seq_len(showCategory)]
    }else{
      if (sum(showCategory %in% names(geneSets)) == 0){
        cli::cli_abort(
          "The value provided to {.arg showCategory} cannot be found in
                {.field SE_data.fgsea@result}."
        )
      }
      geneSets <- geneSets[showCategory]
    }
    names(geneSets) <- .labelBreak(names(geneSets), breaklineN)
    cnetLayout <- ifelse(circular, igraph::layout.circle, igraph::layout_nicely)
    cnetPlot <- ggtangle::cnetplot(geneSets, foldChange=geneList,
      showCategory=length(geneSets),
      node_label=node_label, layout=cnetLayout, cex_label_category=fontSize,
      hilight_alpha = .3) +
      ggplot2::scale_color_gradient2(
        low="#327eba", mid='white', high="#e06663", midpoint=0)
    #cnetPlot <- enrichplot::cnetplot(
    #    res.gsea, foldChange=geneList, node_label=node_label,
    #    cex_label_category=fontSize, layout=cnetLayout,
    #    showCategory=nrow(res.gsea@result))
    if(color){cnetPlot <- cnetPlot +
      ggplot2::guides(colour=ggplot2::guide_colorbar(title=dotColor))}
    TableEdge_cnetPlot <- res.gsea@result |>
        dplyr::select(from=ID, to=core_enrichment) |>
        tidyr::separate_rows(to, sep="/")
    if(color){
        TableNode_cnetPlot <- cnetPlot$data |>
            dplyr::select(id=name, size, foldChange)
    }else{
        TableNode_cnetPlot <- cnetPlot$data |> dplyr::select(id=name, size)
    }
    return(list(cnetPlot=cnetPlot, TableEdge_cnetPlot=TableEdge_cnetPlot,
        TableNode_cnetPlot=TableNode_cnetPlot))
}
