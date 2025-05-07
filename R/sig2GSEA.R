#' @title Gene Set Enrichment Analysis for Signature Interpretation
#' @description Performs pathway enrichment analysis using correlation statistics from signature-transcriptome associations. This core component of the SigFun package enables systematic identification of biological functions associated with gene signatures through GSEA methodology. The function processes SigCor results to rank genes and performs enrichment against specified pathway databases.
#'
#' @param SE_data.cor A `SummarizedExperiment` object containing correlation results from `sigCor`. Must have:
#' - `cor.df` in metadata: Dataframe with correlation statistics (from `sigCor` output)
#' - `rowData`: Gene annotation mapping (ENSEMBL IDs to symbols)
#' @param ranking.method Character. Gene ranking method:
#' - "stat" (default): Use correlation statistics
#' - "pval": Use p-values from correlation analysis
#' @param pathways.all List. Pathway database containing gene sets and associated genes.
#'   Typically obtained via `fgsea::gmtPathways` from MSigDB files (e.g., msigdb.v2023.1.Hs.symbols.gmt)
#' @param NPROC Integer. Number of processors for parallel computation (default=1 for serial processing)
#'
#' @return Returns an augmented `SummarizedExperiment` object with:
#' - `fgseaRes` in metadata: Dataframe containing GSEA results with columns:
#'   * pathway: Pathway name
#'   * pval: Nominal p-value
#'   * padj: Adjusted p-value (FDR)
#'   * ES: Enrichment score
#'   * NES: Normalized enrichment score
#'   * size: Number of genes in pathway
#'   * leadingEdge: Vector of leading edge genes
#' - `DESeq.ranks` in metadata: Named vector of ranked genes
#' - `cor.df` preserved from input
#'
#' @details This function bridges signature-transcriptome correlations with biological interpretation by:
#' 1. Extracting correlation statistics from SigCor results
#' 2. Creating gene rankings based on specified metric (statistics or p-values)
#' 3. Performing GSEA using the fgsea algorithm
#' 4. Storing results in SE object for downstream visualization
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata
#' @importFrom dplyr select rename left_join mutate filter distinct group_by summarize across where
#' @importFrom tibble deframe
#' @importFrom fgsea fgsea
#' @importFrom stats na.omit
#' @export
#'
#' @examples
#' \dontrun{
#' # Load demo data
#' data("demo_GSE181574")
#'
#' # First run correlation analysis
#' SE_cor <- sigCor(SE_data = SE_GSE181574,
#'                 cor.method = "logit",
#'                 Z.transform = FALSE)
#'
#' # Perform GSEA using correlation statistics
#' SE_gsea <- sig2GSEA(
#'   SE_data.cor = SE_cor,
#'   ranking.method = "stat",
#'   pathways.all = pathways.all,
#'   NPROC = 4
#' )
#'
#' # Access results
#' metadata(SE_gsea)$fgseaRes  # GSEA results
#' metadata(SE_gsea)$DESeq.ranks  # Gene rankings
#'
#' # Note: Warnings about p-value estimation are normal for highly significant pathways
#' }
sig2GSEA <- function(SE_data.cor, ranking.method, pathways.all,
                    NPROC=1){
    mapping <- as.data.frame(SummarizedExperiment::rowData(SE_data.cor))
    if ("ensg_id" %in% colnames(mapping)) {
        CodingGene <- mapping %>% dplyr::select(ENSG="ensg_id",
                                            gene_symbol="gene_symbol")
    }else{
        CodingGene <- mapping %>% dplyr::select(ENSG="V2", gene_symbol="V1")
    }
    #
    input <- S4Vectors::metadata(SE_data.cor)$cor.df %>%
        dplyr::select(id=gene, stat=cor)
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
    #data.table::fwrite(DESeq.ranksDF, col.names=TRUE, row.names=FALSE,
    #                    sep='\t', quote=FALSE, na=NA,
    #                    file=file.path(output_path,
    #paste("ranks.notdeframe",gsub("stat","",ranking.method),sep="",".txt")))

    ## run SEA ####
    fgseaRes <- fgsea::fgsea(pathways=pathways.all, stats=DESeq.ranks,
                        maxSize=500L ,minSize=3L ,nproc=NPROC)

    S4Vectors::metadata(SE_data.cor) <- list(fgseaRes=fgseaRes %>%
        dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 5))) %>%
        dplyr::mutate(pathway=gsub('_', " ", x=pathway)),
        cor.df=S4Vectors::metadata(SE_data.cor)$cor.df,
        DESeq.ranks=DESeq.ranks)

    return(SE_data.cor)
}


