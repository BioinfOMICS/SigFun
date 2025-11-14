#' @title ridgePlot
#' @description The \code{ridgePlot} function generates a ridge density
#'   visualization for Gene Set Enrichment Analysis (GSEA) results stored in a
#'   \code{SummarizedExperiment} object. Each ridge represents the distribution
#'   of ranking metrics for genes within a pathway, allowing intuitive
#'   comparison of enrichment signal patterns across multiple gene sets.
#'
#' @param seDataFgsea A \code{SummarizedExperiment} object containing GSEA
#'   results, typically extracted using \code{.extractDF(type = "gseaRaw")}.
#'   The metadata must include enrichment statistics and ranked gene lists.
#' @param showCategory Numeric or character. Specifies the number of top
#'   pathways to display or a vector of specific pathway names to include.
#'   Default is \code{10}.
#' @param breaklineN Integer. Maximum number of characters before inserting a
#'   line break in pathway names for readability. Default is \code{30}.
#' @param fill Character. Variable used for ridge fill color encoding. Must be
#'   one of \code{"pvalue"}, \code{"p.adjust"}, or \code{"qvalue"}.
#'   Default is \code{"p.adjust"}.
#' @param fontSize Numeric. Base font size for axis text and titles.
#'   Default is \code{10}.
#' @param decreasing Logical. Controls the ordering of pathways along the y-axis
#'   based on the median ranking metric. If \code{TRUE}, sorts in descending
#'   order. Default is \code{TRUE}.
#'
#' @details
#'   This function provides a ridge plot visualization inspired by
#'   \pkg{enrichplot::ridgeplot()}, allowing flexible customization and unified
#'   color encoding.
#'
#'   The fill color of each ridge corresponds to -log10(significance) of the
#'   specified \code{fill} variable. Pathways are ordered by the median of their
#'   gene-level ranking metrics, and each curve represents the density of ranked
#'   genes for a given pathway.
#'
#'   Visual encoding:
#'   \itemize{
#'     \item \strong{X-axis:} Gene ranking metric
#'     \item \strong{Y-axis:} Pathway names (auto-wrapped by \code{breaklineN})
#'     \item \strong{Fill color:} -log10(p-value or adjusted value)
#'     \item \strong{Ridge shape:} Distribution of enrichment scores per gene set
#'   }
#'
#' @return A named list with two components:
#'   \describe{
#'     \item{\code{ridgePlot}}{A \code{ggplot2} object displaying ridge density
#'       curves of gene ranking metrics per pathway.}
#'     \item{\code{tableRidgePlot}}{A \code{data.frame} summarizing data used in
#'       the plot, including:
#'       \itemize{
#'         \item \code{original_name}: Original pathway name
#'         \item \code{label_name}: Formatted pathway name with line breaks
#'         \item \code{raw_fill}: Raw significance value (from the \code{fill} column)
#'         \item \code{neg_log10_fill}: -log10 transformed significance
#'         \item \code{ranking_metric}: Gene-level ranking values
#'         \item \code{median_ranking_metric}: Median ranking value per pathway (calculated per group)
#'       }}
#'   }
#'
#' @importFrom DOSE geneInCategory
#' @importFrom ggridges geom_density_ridges
#' @importFrom ggplot2 ggplot aes scale_fill_gradient xlab ylab theme_classic
#'   element_text element_blank element_line scale_y_discrete guide_colorbar
#' @importFrom dplyr mutate select if_else group_by ungroup
#' @importFrom tibble remove_rownames
#' @importFrom rlang check_installed .data
#' @importFrom scales label_number
#'
#' @export
#' @examples
#' data("sig2Fun_result")
#' ridgePlot(seDataFgsea = sig2Fun_result)s
ridgePlot <- function(seDataFgsea, showCategory = 10, breaklineN = 30,
                      fill = "p.adjust", fontSize = 10,
                      decreasing = TRUE) {

    match.arg(fill, c("pvalue", "p.adjust", "qvalue"))
    x <- .extractDF(seDataFgsea, type = 'gseaRaw')

    if (!is(x, "gseaResult"))
        stop("currently only support gseaResult")

    if (!fill %in% colnames(x@result)) {
        stop("'fill' variable not available ...")
    }

    if (inherits(showCategory, 'numeric')) {
        selected <- seq_len(showCategory)
    } else if (inherits(showCategory, "character")) {
        ii <- match(showCategory, x@result$Description)
        if (all(is.na(ii))) {
            ii <- match(showCategory, x@result$ID)
        }
        ii <- ii[!is.na(ii)]
        selected <- x@result[ii, "ID"]
    } else {
        warning("showCategory should be a number of pathways or a vector of selected pathways")
    }

    gs2id <- DOSE::geneInCategory(x)[selected]

    if (x@readable && length(x@gene2Symbol) > 0) {
        id <- match(names(x@geneList), names(x@gene2Symbol))
        names(x@geneList) <- x@gene2Symbol[id]
    }

    gs2val <- lapply(gs2id, function(id) {
        res <- x@geneList[id]
        res <- res[!is.na(res)]
    })

    nn <- names(gs2val)
    i <- match(nn, x$ID)
    nn <- x$Description[i]

    # Remove missing or empty category names to prevent errors in .labelBreak()
    valid_idx <- !is.na(nn)
    nn <- nn[valid_idx]
    gs2val <- gs2val[valid_idx]
    i <- i[valid_idx]

    # Sort pathways by the median of their gene ranking metric
    median_vals <- sapply(gs2val, median, na.rm = TRUE)
    j <- order(median_vals, decreasing = !decreasing)
    len <- sapply(gs2val, length)
    if (length(gs2val) == 0 || all(sapply(gs2val, length) == 0)) {
        cli::cli_abort("No valid gene sets found. Try setting {.arg coreEnrichment = TRUE}.")
    }
    gs2val.df <- data.frame(
        category = rep(nn, times = len),
        color = rep(x[i, fill], times = len),
        value = unlist(gs2val)
    )
    if (!"value" %in% colnames(gs2val.df) || all(is.na(gs2val.df$value))) {
        cli::cli_abort("`ridgePlot()` failed to extract valid ranking metrics. Check input or coreEnrichment parameter.")
    }
    colnames(gs2val.df)[2] <- fill
    gs2val.df$category <- factor(gs2val.df$category, levels = nn[j])

    # plotting
    rlang::check_installed('ggridges', 'for `ridgeplot()`.')
    ridgePlot <- ggplot(gs2val.df,
                        aes(x = .data$value, y = .data$category, fill = -log10(.data[[fill]])) ) +
        ggridges::geom_density_ridges(linewidth = 0.4) +
        scale_fill_gradient(
            low = "#fcae91", high = "#b3200a",
            guide = guide_colorbar(barwidth = 0.9, barheight = 4.5),
            labels = scales::label_number(accuracy = 0.1)
        ) +
        xlab("Gene ranking metric") +
        ylab(NULL) +
        theme_classic(base_size = fontSize) %+replace%
        theme(
            axis.title.x = element_text(
                colour = "black", size = fontSize + 2, face = "bold",
                margin = margin(t = 10)
            ),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.grid.major = element_line(color = "grey95", linewidth = 0.3),
        ) +
        ggplot2::scale_y_discrete(labels = function(x) {
            .labelBreak(x, breaklineN)
        })

    tableRidgePlot <- ridgePlot$data |>
        dplyr::group_by(category) |>
        dplyr::mutate(
            label_name = .labelBreak(category, breaklineN),
            raw_fill = .data[[fill]],
            neg_log10_fill = dplyr::if_else(
                !is.na(.data[[fill]]), -log10(.data[[fill]]), NA_real_
            ),
            median_ranking_metric = median(value, na.rm = TRUE),
            ranking_metric = value
        ) |>
        dplyr::ungroup() |>
        dplyr::select(
            original_name = category, label_name,
            raw_fill, neg_log10_fill,
            ranking_metric, median_ranking_metric
        ) |>
        tibble::remove_rownames()

    return(list(
        ridgePlot = ridgePlot,
        tableRidgePlot = tableRidgePlot
    ))
}
