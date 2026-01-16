#' @title chordPlot
#' @description The \code{chordPlot} function creates a chord diagram
#'   visualization to display the relationships between genes and their
#'   associated pathways from gene set enrichment analysis (GSEA) results. The
#'   chord plot provides an intuitive way to visualize gene-pathway associations,
#'   where genes and pathways are represented as nodes around a circle, connected
#'   by arcs that represent their relationships.
#'
#' @param seDataFgsea A \code{SummarizedExperiment} object containing GSEA
#'   results and pairwise pathway similarity information, typically extracted
#'   using \code{.extractDF(type = "gseaReadable")}.
#' @param showCategory Numeric or character. The number of top pathways to
#'   display, or a vector of specific pathway names to include. Default is
#'   \code{3}.
#' @param breaklineN An integer specifying the maximum number of characters per
#'   line before inserting a line break in pathway names for better readability.
#'   Default is \code{30}.
#' @param fontSize A numeric value between 0 and 1 controlling the font size in
#'   the chord plot. Values closer to 1 result in larger text. Default is
#'   \code{0.5}.
#'
#' @details This function internally calls the helper function \code{.chrod()}
#'   (stored in \code{utils-plot.R}), which uses the \pkg{circlize} package
#'   to create a publication-quality chord diagram. Pathway names are displayed
#'   on the inner side of the circle, while gene names are positioned along the
#'   outer circle, improving readability and visual balance.
#'
#'   The color palette is automatically generated using
#'   \code{viridisLite::viridis()} to ensure a perceptually uniform and
#'   color-blind–friendly representation.
#'
#' @return A named list containing two elements:
#'   \describe{
#'     \item{\code{chordPlot}}{A chord diagram visualization object showing
#'       gene-pathway relationships from the GSEA results. Genes and
#'       pathways are represented as nodes around a circle, connected
#'       by arcs that represent their associations}
#'     \item{\code{tableChordPlot}}{A data.frame containing the chord plot
#'       data with columns: \code{original_name} (original pathway IDs),
#'       \code{label_name} (formatted pathway names with line breaks),
#'       and \code{Gene} (genes associated with each pathway, concatenated
#'       with '/' separator)}
#'   }
#'
#' @importFrom dplyr mutate select group_by distinct
#' @importFrom grDevices dev.cur recordPlot dev.off
#' @importFrom cli cli_abort
#'
#' @export
#' @examples
#' data("sig2Fun_result")
#' chordPlot(seDataFgsea = sig2Fun_result)
#'
#' @note
#' Requires \pkg{circlize} (>= 0.4.15) for chord diagram generation.

# Note: Requires circlize (>= 0.4.15)
chordPlot <- function(
        seDataFgsea, showCategory = 3, breaklineN = 30, fontSize = 0.5) {
    if (isFALSE(is.numeric(fontSize))) {
        cli::cli_abort("The input {.arg fontSize} must be a {.cls numeric}.")
    } else {
        if (fontSize <= 0 | fontSize > 1) {
            cli::cli_abort(
                "The input {.arg fontSize} must be a numeric value between 0
                and 1.")
        }
    }
    data <- .extractDF(seDataFgsea, type = "gseaReadable") |>
        .extractGeneSetsDF(showCategory) |>
        dplyr::mutate(name = .labelBreak(categoryID, breaklineN))
    initial_dev <- grDevices::dev.cur()
    data |> dplyr::select(Gene, name) |> .chrod(fontSize)
    chordPlot <- grDevices::recordPlot()
    if (grDevices::dev.cur() > initial_dev && grDevices::dev.cur() > 1) {
        try(grDevices::dev.off(), silent = TRUE)
    }
    tableChordPlot <- data |> dplyr::group_by(categoryID) |>
        dplyr::mutate(Gene = paste0(Gene, collapse = '/')) |>
        dplyr::select(original_name = categoryID, label_name = name, Gene) |>
        dplyr::distinct()
    return(list(chordPlot = chordPlot, tableChordPlot = tableChordPlot))
}
