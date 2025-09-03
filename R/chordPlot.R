#' @title chordPlot
#' @description The \code{chordPlot} function creates a chord diagram
#' visualization to display the relationships between genes and their
#' associated pathways from gene set enrichment analysis (GSEA) results. The
#' chord plot provides an intuitive way to visualize gene-pathway associations,
#' where genes and pathways are represented as nodes around a circle, connected
#' by arcs that represent their relationships.
#' @param SE_data.fgsea A \code{SummarizedExperiment} object containing GSEA
#'   results in its metadata. The object must contain:
#'   \itemize{
#'     \item \code{metadata(SE_data.fgsea)$gseaResult}: GSEA result object with
#'           pathway information, gene sets, and enrichment statistics
#'   }
#' @param breaklineN An integer specifying the maximum number of characters per
#'   line before inserting a line break in pathway names for better readability.
#'   Default is \code{30}.
#' @param fontSize A numeric value between 0 and 1 controlling the font size in
#' the chord plot. Values closer to 1 result in larger text. Default is
#' \code{0.5}.
#' @param showCategory Either a numeric value specifying the number of top
#'   pathways to display, or a character vector of specific pathway names to
#'   include in the visualization. Default is \code{5}.
#' @return A list containing two elements:
#'   \describe{
#'     \item{\code{chordPlot}}{A chord diagram visualization object showing
#'           gene-pathway relationships from the GSEA results. Genes and
#'           pathways are represented as nodes around a circle, connected
#'           by arcs that represent their associations}
#'     \item{\code{Table_chordPlot}}{A data.frame containing the chord plot
#'           data with columns: \code{original_name} (original pathway IDs),
#'           \code{label_name} (formatted pathway names with line breaks),
#'           and \code{Gene} (genes associated with each pathway, concatenated
#'           with '/' separator)}
#'   }
#' @export
#' @examples
#' data("demo_GSE181574")
#' chordPlot(SE_data.fgsea=GSE181574.sigfun)
chordPlot <- function(
        SE_data.fgsea, showCategory=5, breaklineN=30, fontSize=0.5){
    if(isFALSE(is.numeric(fontSize))){
        cli::cli_abort("The input {.arg fontSize} must be a {.cls numeric}.")
    }else{
        if(fontSize <= 0 | fontSize > 1){
            cli::cli_abort(
                "The input {.arg fontSize} must be a numeric value between 0
                and 1.")
        }
    }
    data <- .extractDF(SE_data.fgsea, type="gseaReadable") |>
        .extractGeneSets(showCategory) |>
        dplyr::mutate(name=.labelBreak(categoryID, breaklineN))
    chordPlot <- data |> dplyr::select(Gene, name) |> .chrod(fontSize)
    Table_chordPlot <- data |> dplyr::group_by(categoryID) |>
        dplyr::mutate(Gene=paste0(Gene, collapse='/')) |>
        dplyr::select(original_name=categoryID, label_name=name, Gene) |>
        dplyr::distinct()
    return(list(chordPlot=chordPlot, Table_chordPlot=Table_chordPlot))
}
