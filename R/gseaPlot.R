#' @title gseaPlot
#' @description The \code{gseaPlot} function creates a GSEA (Gene Set Enrichment
#' Analysis) plot visualization from GSEA results data. This function serves as
#' a wrapper around \code{enrichplot::gseaplot2} to simplify the process of
#' generating enrichment score plots. The plot displays the running enrichment
#' score, gene positions, and statistical information for selected gene sets,
#' providing a comprehensive view of how genes within a pathway contribute to
#' the enrichment signal.
#' @param SE_data.fgsea A \code{SummarizedExperiment} object or data structure
#'   containing GSEA results from fgsea analysis. The object must contain:
#'   \itemize{
#'     \item GSEA result data that can be processed by the internal
#'           \code{.extractDF} function with type "gseaRaw"
#'     \item Enrichment statistics including normalized enrichment scores,
#'           p-values, and gene set information
#'   }
#' @param showCategory Either a numeric value specifying the number of
#' top-ranked gene sets to display (based on statistical significance), or a
#' character vector specifying the exact names of gene sets to visualize. When
#' numeric, the function selects the first N gene sets from the results.
#' Default is \code{10}.
#' @return A ggplot2 object containing the GSEA plot visualization showing the
#' running enrichment score curves, gene positions along the ranked gene list,
#' and statistical information for the selected gene sets. The plot provides a
#' comprehensive view of how genes within each pathway contribute to the
#' enrichment signal.
#' @export
#' @examples
#' data("demo_GSE181574")
#' gseaPlot(SE_data.fgsea=GSE181574.sigfun)
gseaPlot <- function(SE_data.fgsea, showCategory=5){
    res.gsea <- .extractDF(SE_data.fgsea, type="gseaRaw")
    if (is.numeric(showCategory)) {
        geneSetID <- seq_len(showCategory)
    }else{
        geneSetID <- showCategory
    }
    gseaPlot <- enrichplot::gseaplot2(res.gsea, geneSetID=geneSetID)
    return(gseaPlot)
}
