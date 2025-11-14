#' @title Signature-to-Function Analysis
#' @description Analyzes gene signatures and identifies their associated
#' biological functions by using the whole transcriptome as a surrogate.
#' This function provides a comprehensive solution for functional analysis of
#' multi-gene signatures regardless of their size, addressing limitations of
#' conventional methods like ORA and GSEA for small gene sets. It can handle
#' signatures in various formats including binary values (e.g., high/low risk
#' classifications) and continuous values (e.g., prognostic risk scores).
#'
#' @param seData A SummarizedExperiment object containing:
#'        - assays: Gene expression matrix, genes in rows, samples in columns
#'        - colData: Signature scores for each sample
#'        - rowData: Gene annotation including gene symbols that match the
#'        pathway database
#'        See vignette for detailed structure requirements.
#' @param rankingMethod Character. Method to rank genes for pathway analysis:
#'        - "cor" (default): Use test statistics
#'        - "pval": Use p-values
#' @param species Character. Organism of the dataset:
#'        - "human" (default)
#'        - "mouse"
#' @param corMethod Character. Statistical method for correlation analysis:
#'        - "spearman" (default): For general correlations
#'        - "pearson": For linear relationships
#'        - "kendall": For non-parametric rank correlation
#'        - "logit": For binary signatures (e.g., high/low risk classification)
#' @param t2g List containing biological functions and their
#' corresponding genes.
#'        Should be obtained from MSigDB (e.g., msigdb.v2023.1.Hs.symbols.gmt)
#'        or a custom gene set.
#'        Each list element represents a pathway with genes as vector elements.
#' @param topN Integer. Number of top significant functions to display in output
#' plots.
#'        Default is 10.
#' @param zTransform Logical. Whether to perform z-transformation on gene
#' expression data:
#'        - FALSE (default): No transformation
#'        - TRUE: Apply z-transformation (standardization)
#'        For binary signatures (e.g., using logit correlation), set to FALSE.
#' @param significantType Character. Method for statistical significance
#' filtering:
#'        - "pval": Use p-values (default)
#'        - "qval": Use Q-values (adjusted p-values)
#' @param strings Character vector. Specifies which pathway categories to
#' include in visualization.
#'        Default includes: c("GOBP", "GOCC", "GOMF", "KEGG", "REACTOME",
#'        "WP", "HALLMARK", "SIGNALING")
#'        For example, "HALLMARK" will filter for pathways containing it.
#'
#' @return Returns an enhanced SummarizedExperiment object with additional slots
#' containing:
#'         - Correlation results between signature and gene expression
#'         - Pathway enrichment analysis results
#'         - Visualization objects:
#'           * barplots: Bar plots showing top functions associated with the
#'           signature, ordered by significance and colored by association
#'           direction
#'           * heatmap: Heatmaps displaying NES values, p-values, score
#'           distribution, and names of significantly associated functions
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom stats binomial glm cor
#' @importFrom fgsea gmtPathways
#' @importFrom methods new
#' @importFrom grDevices pdf dev.off
#' @importFrom utils write.csv
#' @include utility.R
#' @include sigCor.R
#' @include plot_heat.R
#' @export
#'
#' @examples
#' # Load demo dataset
#' data("expr.data")
#' data("mapping")
#' data("SIG_MAT")
#' data("t2g")
#' seData <- SummarizedExperiment::SummarizedExperiment(
#' assays=list(abundance=as.matrix(expr.data)),
#' rowData=S4Vectors::DataFrame(mapping, row.names=mapping$ensg_id),
#' colData=S4Vectors::DataFrame(SIG_MAT))
#' # For analysis of binary signature (e.g., MammaPrint high/low risk)
#' res <- sig2Fun(seData=seData, t2g=t2g, rankingMethod="cor",
#' species="human", corMethod="logit", topN=10, zTransform=FALSE,
#' significantType="pval", strings=c("HALLMARK"))
#' # To visualize Heatmap results: res$heatmap
sig2Fun <- function(seData, t2g, rankingMethod="cor", topN=10,
                    species="human", corMethod="spearman", zTransform=FALSE,
                    significantType="pval",
                    strings=c("GOBP","REACTOME","HALLMARK","SIGNALING")) {
    .classCheck(seData, "SummarizedExperiment")
    .classCheck(t2g, "data.frame")
    #pathway.name <- unique(t2g$gs_name)
    #new_pathways.all <- lapply(pathway.name, function(x){
    #    res <- t2g %>% dplyr::filter(gs_name==x) %>%
    #    dplyr::select(ensembl_gene) %>%
    #    dplyr::pull()
    #    return(res)
    #})
    #names(new_pathways.all) <- pathway.name
    new_pathways.all <- split(t2g$ensembl_gene, t2g$gs_name)
    seDataCor <- seData
    if(!("cor.df" %in% names(S4Vectors::metadata(seDataCor)))){
        seDataCor <- sigCor(seData=seData, corMethod=corMethod,
                            zTransform=zTransform)
        S4Vectors::metadata(seDataCor) <-
            list(cor.df=S4Vectors::metadata(seDataCor)$cor.df)
    }

    geneList <- seDataCor@metadata$cor.df %>%
        dplyr::select(dplyr::all_of(rankingMethod)) %>% dplyr::pull()
    names(geneList) <- seDataCor@metadata$cor.df$gene
    geneList <- sort(geneList, decreasing=TRUE)
    clusterProfiler <- clusterProfiler::GSEA(geneList, TERM2GENE=t2g)

    seDataFgsea <- seDataCor
    S4Vectors::metadata(seDataFgsea) <- list(gseaResult=clusterProfiler,
                                             cor.df=S4Vectors::metadata(seDataCor)$cor.df)

    # heatmap
    heatmap.list <- NULL
    for (i in seq_along(strings)) {
        heatmap.list[[strings[i]]] <- plot_heat(
            seDataFgsea = seDataFgsea,
            pathwaysAll = new_pathways.all,
            significantType = significantType,
            strings = strings[i],
            topN = topN,
            rankingMethod = rankingMethod
        )
    }
    seDataFgsea@metadata$heatmap  <- heatmap.list

    return(seDataFgsea)
}
