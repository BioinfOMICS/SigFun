#' @title Visualization of Pathway Enrichment Results
#' @description Generates heatmap visualizations of Gene Set Enrichment Analysis (GSEA) results,
#' displaying normalized enrichment scores (NES), significance values, and gene score distributions.
#' This function creates a multi-panel heatmap showing:
#' \itemize{
#'   \item{NES values for significant pathways}
#'   \item{Statistical significance (-log10 transformed p/q-values)}
#'   \item{Gene score distribution patterns}
#'   \item{Pathway names with leading edge genes}
#' }
#' Color gradients indicate association direction (red = positive, blue = negative).
#'
#' @param SE_data.fgsea A `SummarizedExperiment` object containing:
#' \itemize{
#'   \item{`cor.df` in metadata: Correlation statistics from `sigCor`}
#'   \item{`fgseaRes` in metadata: GSEA results from `sig2GSEA`}
#' }
#' @param significat_type Character. Significance metric for filtering:
#' \itemize{
#'   \item{"pval" (default): Use nominal p-values}
#'   \item{"qval": Use false discovery rate adjusted q-values}
#' }
#' @param topN Integer. Number of top significant pathways to display per
#' category (default=10).
#' @param pathways.all List. Pathway database from MSigDB
#' (e.g. msigdb.v2023.1.Hs.symbols.gmt).
#' @param strings Character vector. Pathway categories to visualize
#' (e.g. "HALLMARK", "KEGG").
#' Default includes major MSigDB collections.
#' @param ranking.method Character. Ranking metric used in GSEA:
#' \itemize{
#'   \item{"core" (default): Correlation statistics}
#'   \item{"pval": P-values from correlation analysis}
#' }
#' @param colwidths Vector of five elements corresponding to column width for
#' grid.arrange
#' @param breaklineN An integer specifying the number of characters after which
#' to break pathway names into multiple lines for better readability.
#' Default is \code{30}.
#' @param fontSize Numeric value used to control pathway label font size.
#' @return Returns a list of ggplot objects containing heatmaps, with one
#' element per category specified in `strings`. Each heatmap contains four
#' panels:
#' \enumerate{
#'   \item{NES values}
#'   \item{Significance values (-log10 transformed)}
#'   \item{Score distribution patterns}
#'   \item{Pathway names with leading edge genes}
#' }
#'
#' @examples
#' \dontrun{
#' # After running sig2Fun analysis
#' data("demo_GSE181574")
#' res <- sig2Fun(SE_GSE181574, strings = "HALLMARK")
#'
#' # Generate heatmap for HALLMARK pathways
#' hm <- plot_heat(
#'   SE_data.fgsea = res,
#'   strings = "HALLMARK",
#'   pathways.all = pathways.all,
#'   topN = 5
#' )
#'
#' # Access heatmap plot
#' hm$HALLMARK
#' }
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata
#' @importFrom dplyr filter select left_join mutate arrange desc slice group_by summarize distinct
#' @importFrom ggplot2 ggsave
#' @export

plot_heat <- function(
    SE_data.fgsea, significant_type="pval", topN=10, pathways.all,
    strings=c("GOBP", "GOCC", "GOMF", "KEGG", "REACTOME", "WP", "HALLMARK"),
    ranking.method="cor", colwidths=c(.8, .8, .8, 3.6, 5),
    breaklineN=30, fontSize=8) {
    if(length(colwidths) != 5){
      cli::cli_abort(
        'The {.arg colwidths} must be a vector of length 5, which controls the
        width of "NES bar", "NES", "P-value", "Enrichment plot", and
        "Pathway name".'
      )
    }
    cor.df <- S4Vectors::metadata(SE_data.fgsea)$cor.df
    RES_GSEA <- S4Vectors::metadata(SE_data.fgsea)$gseaResult %>%
      tibble::as_tibble() %>%
        dplyr::slice(grep(paste(strings, collapse="|"), ID)) %>%
        dplyr::mutate(`-log10(pvalue)`=-log10(pvalue))
    if (significant_type == "pval") {
        RES_GSEA <- RES_GSEA %>% dplyr::filter(pvalue < 0.05)
    } else if (significant_type == "padj") {
        RES_GSEA <- RES_GSEA %>% dplyr::filter(p.adjust < 0.05)
    }
    CodingGene <- SummarizedExperiment::rowData(SE_data.fgsea) %>%
        as.data.frame() %>% dplyr::filter(gene_biotype == "protein_coding") %>%
        dplyr::select(ENSG=ensg_id, gene_symbol=gene_symbol)
    ranksDF <- cor.df %>% dplyr::filter(!is.na(gene)) %>%
        dplyr::select(ENSG=gene, stat=ranking.method) %>%
        dplyr::left_join(CodingGene, by="ENSG") %>%
        dplyr::mutate(abs.stat=abs(stat)) %>%
        dplyr::select(ENSG, stat) %>%
        na.omit() %>% dplyr::distinct() %>% dplyr::group_by(ENSG) %>%
        dplyr::summarize(ranking.method=mean(as.numeric(stat))) %>%
        dplyr::filter(ranking.method != 'Inf' & ranking.method != '-Inf') %>%
        tibble::deframe()
    tmp.res <- lapply(strings, function(str) {
        RES_NES <- RES_GSEA %>% dplyr::slice(grep(str, ID))
        RES_NES_top10 <- RES_NES %>%
            dplyr::arrange(dplyr::desc(NES)) %>%
            dplyr::filter(NES >= 0)
        if(nrow(RES_NES_top10) >= topN){
            RES_NES_top10 <- RES_NES_top10 %>%
                dplyr::slice(seq_len(topN))
        }
        RES_NES_bottom10 <- RES_NES %>%
            dplyr::arrange(dplyr::desc(-NES)) %>%
            dplyr::filter(NES <= 0)
        if(nrow(RES_NES_bottom10) >= topN){
            RES_NES_bottom10 <- RES_NES_bottom10 %>%
                dplyr::slice(seq_len(topN))
        }
        RES_NES <- rbind(RES_NES_top10, RES_NES_bottom10) %>%
            dplyr::arrange(dplyr::desc(NES))

        tmp <- .plotGseaTable(
          pathways=pathways.all[RES_NES$ID], stats=ranksDF, fgseaRes=RES_GSEA,
          colwidths=colwidths, breaklineN=breaklineN, fontSize=fontSize)
        tmp
    })
    names(tmp.res) <- strings
    return(tmp.res)
}



