#' @title plot_heat
#' @description  This is a plot component, suggesting the package performs Gene
#' Set Enrichment Analysis (GSEA). It uses the SigCor results (correlation
#' statistics) and/or p-values, along with the ontology profile, to identify
#' enriched pathways or gene sets associated with the input signature.
#' @param SE_data.fgsea Summary Experiment Object. In the metadata of the input,
#' it must be a Summary Experiment Object and should contain the following two
#' objects "cor.df" and "fgseaRes". "cor.df" is the result of SigCor. Add the
#' cor.df to the metadata by running SigCor. "fgseaRes" is the result of
#' Sig2GSEA. Add the fgseaRes to the metadata by running Sig2GSEA
#' @param output_path Character. Location of the output path.
#' @param significat_type Character. Filter or statistical significance.
#' "pval" for p-value and "qval" for Q-value
#' @param strings vector. Pathways header want to be plotted.
#' @param pathways.all gmt file. loading by fgsea::gmtPathways. Should be
#' obtained from MSigDB: msigdb.v2023.1.Hs.symbols.gmt
#' @param ranking.method Character. "stat" (default) or "pval"
#' @param plot_out logic. Plot the heatmap to output_path (Default: FALSE).
#' @param topN Integer. Number of significant functions in hte output bar plots.
#' Default is 10.
#' @return Return heatmaps to the outpath specified.
#' @export
#' @examples
#' data(GSEA_data)
#' output_path <- tempdir()
#' heatmap <- plot_heat(SE_data.fgsea = GSEA_data, output_path = output_path,
#' pathways.all = pathways.all, significat_type = "pval", strings = c("KEGG"),
#' ranking.method = "stat", plot_out=FALSE,topN=3)

plot_heat <- function(SE_data.fgsea, output_path, significat_type="pval",
                strings=c("GOBP","GOCC","GOMF","KEGG","REACTOME","WP"),
                topN=10,pathways.all, ranking.method="stat",plot_out=TRUE){
    cor.df <- SE_data.fgsea@metadata$cor.df
    RES_GSEA <- SE_data.fgsea@metadata$fgseaRes %>% dplyr::filter(pval<0.05)
    RES_NES_total <- RES_GSEA %>% dplyr::slice(
        grep(paste0("[",paste(paste(strings, " "), collapse = "|"),"]"),
            pathway)) %>% dplyr::mutate(`-log10(pvalue)` = -log10(pval))
    RES_NES_total$leadingEdge <- gsub('[|]',';',RES_NES_total$leadingEdge)
    RES_NES_total$pathway <- gsub(" ","_",RES_NES_total$pathway)
    mapping <- as.data.frame(SummarizedExperiment::colData(SE_data.fgsea))
    mapping.tmp <- mapping %>% dplyr::select(gene=gene_symbol, ensg=ensg_id)
    res.out.ensg <- dplyr::left_join(cor.df %>% dplyr::select(ensg=gene,cor),
                                    mapping.tmp) %>% dplyr::filter(!is.na(ensg))
    CodingGene <- mapping %>% dplyr::filter(gene_biotype=="protein_coding") %>%
        dplyr::select(ENSG=ensg_id,gene_symbol=gene_symbol)
    ranksDF  <-  res.out.ensg %>% dplyr::select(id=ensg,stat=cor)  %>%
        dplyr::rename(ENSG=id) %>% dplyr::select(ENSG,stat) %>%
        dplyr::left_join(CodingGene, by ="ENSG") %>%
        dplyr::mutate(abs.stat=abs(stat))%>%
        dplyr::select(gene_symbol, all_of(ranking.method)) %>%
        na.omit() %>% dplyr::distinct() %>% dplyr::group_by(gene_symbol) %>%
    dplyr::summarize(ranking.method=mean(as.numeric(get(ranking.method)))) %>%
        dplyr::filter(ranking.method!='Inf'&ranking.method!='-Inf') %>%
        tibble::deframe()
    tmp.res <- NULL
    for (j in seq_along(strings)) {
        RES_NES_top10 <- RES_NES_total %>%
            dplyr::slice(grep(paste0(strings[j], "_"),RES_NES_total$pathway))%>%
            dplyr::filter(.data[[significat_type]]<0.05) %>%
            dplyr::arrange(dplyr::desc(NES)) %>% dplyr::slice(seq_len(topN))
        RES_NES_bottom10 <- RES_NES_total %>%
        dplyr::slice(grep(paste0(strings[j], "_"), RES_NES_total$pathway)) %>%
        dplyr::filter(.data[[significat_type]]<0.05) %>%
        dplyr::arrange(dplyr::desc(-NES)) %>%
        dplyr::slice(seq_len(topN)) %>% dplyr::arrange(dplyr::desc(NES))
        tmp <- .plotGseaTable(pathways.all, Pathways=pathways.all[
            c(RES_NES_top10$pathway,RES_NES_bottom10$pathway)],
            stats = ranksDF, fgseaRes = RES_NES_total)
        if(plot_out){ggplot2::ggsave(filename = file.path(output_path,
        paste0("GSEA_heatmap_",strings[j],".png")),
        plot = tmp, width = 15, height = 10)}
        tmp.res[[strings[j]]] <- tmp
    }
    return(tmp.res)
}



