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
#' @param strings Character vector. Pathway categories to visualize (e.g. "HALLMARK", "KEGG").
#' Default includes major MSigDB collections.
#' @param topN Integer. Number of top significant pathways to display per category (default=10).
#' @param pathways.all List. Pathway database from MSigDB (e.g. msigdb.v2023.1.Hs.symbols.gmt).
#' @param ranking.method Character. Ranking metric used in GSEA:
#' \itemize{
#'   \item{"stat" (default): Correlation statistics}
#'   \item{"pval": P-values from correlation analysis}
#' }
#'
#' @return Returns a list of ggplot objects containing heatmaps, with one element per category
#' specified in `strings`. Each heatmap contains four panels:
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

plot_heat <- function(SE_data.fgsea, significat_type="pval",
    strings=c("GOBP", "GOCC", "GOMF", "KEGG", "REACTOME", "WP", "HALLMARK"),
    topN=10, pathways.all, ranking.method="stat", plot_out=TRUE) {
    cor.df <- S4Vectors::metadata(SE_data.fgsea)$cor.df
    RES_GSEA <- S4Vectors::metadata(SE_data.fgsea)$fgsea %>% dplyr::filter(pval < 0.05)
    RES_NES_total <- RES_GSEA %>% dplyr::slice(grep(
    paste0("[", paste(paste(strings, " "), collapse="|"), "]"), pathway)) %>%
    dplyr::mutate(`-log10(pvalue)`=-log10(pval))
    RES_NES_total$leadingEdge <- gsub('[|]', ';', RES_NES_total$leadingEdge)
    RES_NES_total$pathway <- gsub(" ", "_", RES_NES_total$pathway)
    mapping <- as.data.frame(SummarizedExperiment::rowData(SE_data.fgsea))
    mapping.tmp <- mapping %>% dplyr::select(gene=gene_symbol, ensg=ensg_id)
    res.out.ensg <- dplyr::left_join(cor.df %>% dplyr::select(ensg=gene, cor),
        mapping.tmp) %>% dplyr::filter(!is.na(ensg))
    CodingGene <- mapping %>% dplyr::filter(gene_biotype == "protein_coding")%>%
        dplyr::select(ENSG=ensg_id, gene_symbol=gene_symbol)
    ranksDF <- res.out.ensg %>% dplyr::select(id=ensg, stat=cor) %>%
        dplyr::rename(ENSG=id) %>% dplyr::select(ENSG, stat) %>%
        dplyr::left_join(CodingGene, by="ENSG") %>%
        dplyr::mutate(abs.stat=abs(stat)) %>%
        dplyr::select(gene_symbol, all_of(ranking.method)) %>%
        na.omit() %>% dplyr::distinct() %>% dplyr::group_by(gene_symbol) %>%
    dplyr::summarize(ranking.method=mean(as.numeric(get(ranking.method)))) %>%
    dplyr::filter(ranking.method != 'Inf' & ranking.method != '-Inf') %>%
        tibble::deframe()
    tmp.res <- lapply(strings, function(str) {
        RES_NES_top10 <- RES_NES_total %>%
            #dplyr::slice(grep(paste0(str, "_"), RES_NES_total$pathway)) %>%
            dplyr::slice(grep(str, RES_NES_total$pathway)) %>%
            dplyr::filter(.data[[significat_type]] < 0.05) %>%
            dplyr::arrange(dplyr::desc(NES)) %>% dplyr::slice(seq_len(topN))
        RES_NES_bottom10 <- RES_NES_total %>%
    #dplyr::slice(grep(paste0(str, "_"), RES_NES_total$pathway)) %>%
    dplyr::slice(grep(str, RES_NES_total$pathway)) %>%
    dplyr::filter(.data[[significat_type]] < 0.05) %>%
    dplyr::arrange(dplyr::desc(-NES)) %>% dplyr::slice(seq_len(topN)) %>%
        dplyr::arrange(dplyr::desc(NES))
        tmp <- .plotGseaTable(pathways.all,
    Pathways=pathways.all[c(RES_NES_top10$pathway, RES_NES_bottom10$pathway)],
    stats=ranksDF,fgseaRes=RES_NES_total)
        #if (plot_out) {
        #ggplot2::ggsave(plot=tmp,width=15,height=10,
        #filename=file.path(output_path, paste0("GSEA_heatmap_", str, ".png")))
        #}
        tmp
    })
    names(tmp.res) <- strings
    return(tmp.res)
}



