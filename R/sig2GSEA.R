#' @title sig2GSEA
#' @description  This is a key component, suggesting the package performs Gene
#' Set Enrichment Analysis (GSEA). It uses the SigCor results (correlation
#' statistics) and/or p-values, along with the ontology profile, to identify
#' enriched pathways or gene sets associated with the input signature.
#' @param SE_data.cor Summarized Experiment Object. In the metadata of the input
#' , it must be a Summary Experiment Object and should contain the following
#' objects "cor.df". "cor.df" is the result of SigCor.
#' Add the cor.df to the metadata by running SigCor.
#' @param ranking.method Character. "stat" (default) or "pval"
#' @param pathways.all gmt file. loading by fgsea::gmtPathways
#' @param output_path Character. Location of the output path.
#' @param NPROC Parallelization parameter used in bplapply (default=1).
#' Integer Number of permutations to do (default=1000).
#' @return Return a SummarizedExperiment object containing analysis results.
#' SigCor will add "fgseaRes" to the metadata of the input object which is the
#' necessary input for plot_bar and plot_heat.
#' @include utility.R
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr select rename left_join mutate filter distinct group_by
#'   summarize across where
#' @importFrom tibble deframe
#' @importFrom data.table fwrite
#' @importFrom fgsea fgsea
#' @importFrom stats na.omit
#' @export
#' @examples
#' data(sig2GSEA_input)
#' temp_dir <- tempdir()
#' result <- sig2GSEA(SE_data.cor=test_SE, ranking.method="stat",
#' output_path=temp_dir, pathways.all=test_pathways, NPROC=2)

sig2GSEA <- function(SE_data.cor, ranking.method, output_path, pathways.all,
                    NPROC=1){
    mapping <- as.data.frame(SummarizedExperiment::colData(SE_data.cor))
    if ("ensg_id" %in% colnames(mapping)) {
        CodingGene <- mapping %>% dplyr::select(ENSG="ensg_id",
                                            gene_symbol="gene_symbol")
    }else{
        CodingGene <- mapping %>% dplyr::select(ENSG="V2", gene_symbol="V1")
    }
    #
    input <- metadata(SE_data.cor)$cor.df %>% dplyr::select(id=gene, stat=cor)
    input <- input %>% dplyr::rename(ENSG=id) %>% dplyr::select(ENSG,stat) %>%
        dplyr::left_join(CodingGene, by ="ENSG" )

    DESeq.ranksDF <- input%>%
        dplyr::mutate(abs.stat=abs(stat))%>%
        dplyr::select(gene_symbol, all_of(ranking.method)) %>%
        na.omit() %>% dplyr::distinct() %>% dplyr:: group_by(gene_symbol) %>%
        dplyr::summarize(ranking.method=
            mean(as.numeric(get(ranking.method)))) %>%
        dplyr::filter(ranking.method!='Inf'&ranking.method!='-Inf')
    DESeq.ranks <- DESeq.ranksDF %>% tibble::deframe()
    assign('DESeq.ranks',DESeq.ranks)
    data.table::fwrite(DESeq.ranksDF, col.names=TRUE, row.names=FALSE,
                        sep='\t', quote=FALSE, na=NA,
                        file=file.path(output_path,
    paste("ranks.notdeframe",gsub("stat","",ranking.method),sep="",".txt")))

    ## run SEA ####
    fgseaRes <- fgsea::fgsea(pathways=pathways.all, stats=DESeq.ranks,
                        maxSize=500L ,minSize=3L ,nproc=NPROC)

    metadata(SE_data.cor) <- list(fgseaRes=fgseaRes %>%
        dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 5))) %>%
        dplyr::mutate(pathway=gsub('_', " ", x=pathway)))

    return(SE_data.cor)
}


