#' @title Signature-to-Function Analysis
#' @description Analyzes gene signatures and identifies their associated
#' biological functions by using the whole transcriptome as a surrogate.
#' This function provides a comprehensive solution for functional analysis of
#' multi-gene signatures regardless of their size, addressing limitations of
#' conventional methods like ORA and GSEA for small gene sets. It can handle
#' signatures in various formats including binary values (e.g., high/low risk
#' classifications) and continuous values (e.g., prognostic risk scores).
#'
#' @param SE_data A SummarizedExperiment object containing:
#'        - assays: Gene expression matrix, genes in rows, samples in columns
#'        - colData: Signature scores for each sample
#'        - rowData: Gene annotation including gene symbols that match the
#'        pathway database
#'        See vignette for detailed structure requirements.
#' @param ranking.method Character. Method to rank genes for pathway analysis:
#'        - "stat" (default): Use test statistics
#'        - "pval": Use p-values
#' @param species Character. Organism of the dataset:
#'        - "human" (default)
#'        - "mouse"
#' @param cor.method Character. Statistical method for correlation analysis:
#'        - "spearman" (default): For general correlations
#'        - "pearson": For linear relationships
#'        - "kendall": For non-parametric rank correlation
#'        - "logit": For binary signatures (e.g., high/low risk classification)
#' @param pathways_all List containing biological functions and their
#' corresponding genes.
#'        Should be obtained from MSigDB (e.g., msigdb.v2023.1.Hs.symbols.gmt)
#'        or a custom gene set.
#'        Each list element represents a pathway with genes as vector elements.
#' @param topN Integer. Number of top significant functions to display in output
#' plots.
#'        Default is 10.
#' @param Z.transform Logical. Whether to perform z-transformation on gene
#' expression data:
#'        - FALSE (default): No transformation
#'        - TRUE: Apply z-transformation (standardization)
#'        For binary signatures (e.g., using logit correlation), set to FALSE.
#' @param significat_type Character. Method for statistical significance
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
#' @include sig2GSEA.R
#' @include plot_bar.R
#' @include plot_heat.R
#' @export
#'
#' @examples
#' # Load demo dataset
#' data("demo_GSE181574")
#'
#' # For analysis of binary signature (e.g., MammaPrint high/low risk)
#' res <- sig2Fun(
#'   SE_data = SE_GSE181574,
#'   ranking.method = "stat",
#'   species = "human",
#'   cor.method = "logit",  # Use logit for binary variables
#'   pathways_all = pathways.all,
#'   topN = 10,
#'   Z.transform = FALSE,  # No normalization needed for binary variables
#'   significat_type = "pval",
#'   strings = c("HALLMARK")
#' )
#'
#' # To visualize results:
#' # 1. Bar plot: res$barplots
#' # 2. Heatmap: res$heatmap

sig2Fun <- function(SE_data, ranking.method="stat", species="human",
                    cor.method="spearman", pathways_all,
                    topN=10, Z.transform=FALSE, significat_type="pval",
                    strings=c("GOBP","GOCC","GOMF","KEGG","REACTOME","WP","HALLMARK","SIGNALING")) {
    .classCheck(SE_data, "SummarizedExperiment")
    .classCheck(pathways_all, "list")

    SE_data.cor <- SE_data
        if(!("cor.df" %in% names(metadata(SE_data.cor)))){
            SE_data.cor <- sigCor(SE_data=SE_data, cor.method=cor.method,
             Z.transform=Z.transform)
            metadata(SE_data.cor) <- list(cor.df=metadata(SE_data.cor)$cor.df)
        }
    SE_data.fgsea <- SE_data.cor
        if(!("fgseaRes" %in% names(metadata(SE_data.fgsea)))){
            SE_data.fgsea <- sig2GSEA(SE_data.cor=SE_data.cor,
                ranking.method=ranking.method,
                pathways.all=pathways_all)
            metadata(SE_data.fgsea) <- list(fgsea=metadata(SE_data.fgsea)$fgsea,
                                            cor.df=metadata(SE_data.cor)$cor.df)
        }

    barplots <- plot_bar(SE_data.fgsea=SE_data.fgsea,
    topN=topN, significat_type=significat_type, strings=strings)
    SE_data.fgsea$barplots <- barplots

    heatmap <- plot_heat(SE_data.fgsea=SE_data.fgsea,
        strings=strings, significat_type=significat_type, topN=topN,
        pathways.all=pathways_all, ranking.method=ranking.method)
    SE_data.fgsea$heatmap  <- heatmap

    return(SE_data.fgsea)
}
