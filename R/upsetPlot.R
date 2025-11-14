#' @title upsetPlot
#' @description The \code{upsetPlot} function generates an UpSet plot for
#'   visualizing overlapping genes across multiple pathways from gene set
#'   enrichment analysis (GSEA) results. The plot combines two panels: the
#'   upper panel shows either a boxplot or barplot summarizing pathway
#'   combinations, while the lower panel displays an intersection matrix
#'   illustrating pathway overlaps.
#'
#' @param seDataFgsea A \code{SummarizedExperiment} object containing GSEA
#'   results and correlation data, typically extracted using
#'   \code{.extractDF(type = "gseaRaw")}.
#' @param showCategory Numeric or character. Either a number specifying how many
#'   top pathways to display, or a vector of pathway names to include. Default
#'   is \code{5}.
#' @param breaklineN Integer. Number of characters after which to break long
#'   pathway names into multiple lines for readability. Default is \code{30}.
#' @param type Character. Type of plot to display in the upper panel. Must be
#'   one of \code{"bar"} or \code{"box"}. Default is \code{"box"}.
#' @param fontSize Numeric. Text size for axis labels and annotations in both
#'   panels. Default is \code{10}.
#' @param fillColor Character. Color used to fill bars, boxes, or matrix cells
#'   in the plot. Default is \code{"#5979A3"}.
#'
#' @details
#'   This function is designed to visualize the overlap between gene sets from
#'   GSEA results using the UpSet plot framework. It provides a compact
#'   alternative to traditional Venn diagrams, particularly suitable for
#'   comparing multiple pathways.
#'
#'   The upper panel summarizes statistical patterns among pathway combinations:
#'   \itemize{
#'     \item \strong{Boxplot} (\code{type = "box"}): Shows the distribution of
#'       correlation coefficients for genes shared across pathway combinations.
#'       Requires correlation data in the input object.
#'     \item \strong{Barplot} (\code{type = "bar"}): Displays the count of
#'       overlapping genes for each pathway intersection.
#'   }
#'
#'   The lower panel represents an intersection matrix where:
#'   \itemize{
#'     \item Horizontal lines connect dots to show which pathways participate
#'       in each combination
#'     \item Dot positions indicate individual pathway membership
#'     \item The matrix allows quick identification of complex overlaps
#'   }
#'
#' @return A named list containing two components:
#'   \describe{
#'     \item{\code{upsetPlot}}{A combined \code{ggplot2} object (via
#'       \code{ggpubr::ggarrange}) containing two aligned panels: the upper
#'       summary panel (boxplot or barplot) and the lower intersection matrix.}
#'     \item{\code{tableUpsetPlot}}{A \code{data.frame} containing the raw data
#'       used for generating the plot, including pathway combinations, overlap
#'       counts, gene identities, and correlation statistics (when applicable).}
#'   }
#'
#' @importFrom ggpubr ggarrange
#' @importFrom dplyr select filter
#' @importFrom ggplot2 ggplot aes geom_point geom_segment theme_minimal theme
#'   element_blank element_text element_line scale_y_discrete labs
#'
#' @export
#' @examples
#' data("sig2Fun_result")
#' upsetPlot(seDataFgsea = sig2Fun_result)
#'
#' @note
#' Requires \pkg{ggpubr} (>= 0.6.0) and \pkg{dplyr} (>= 1.1.0) for full compatibility.

# Note: Requires ggpubr (>= 0.6.0) and dplyr (>= 1.1.0)
upsetPlot <- function(seDataFgsea, showCategory = 5, breaklineN = 30, type = "box",
                      fontSize = 10, fillColor = "#5979A3") {

    match.arg(type, c("bar", "box"))

    gseaRaw <- .extractDF(seDataFgsea, type = "gseaRaw")
    geneSets <- .extractGeneSets(gseaRaw, showCategory)
    upsetData <- .upsetData(geneSets, gseaRaw, breaklineN, type)

    upsetTop <- switch(
        type,
        bar = .upsetBarplot(upsetData$data, fontSize = fontSize, fillColor = fillColor),
        box = .upsetBoxplot(upsetData$data, fontSize = fontSize, fillColor = fillColor)
    )

    upsetBottom <- .upsetPlot(upsetData$data, fontSize = fontSize, fillColor = fillColor)

    nPath <- length(unique(geneSets$categoryID))
    if (nPath < 5) {
        heights <- c(1 - nPath / 20, nPath / 20)
    } else {
        heights <- c(0.5, 0.5)
    }

    upsetPlot <- ggpubr::ggarrange(
        upsetTop, upsetBottom, ncol = 1, align = "v", heights = heights
    )

    tableUpsetPlot <- upsetData$Rawdata

    return(list(upsetPlot = upsetPlot, tableUpsetPlot = tableUpsetPlot))
}
