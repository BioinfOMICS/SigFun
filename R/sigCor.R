#' @title Signature-Transcriptome Correlation Analysis
#' @description Calculates genome-wide correlations between a clinical signature
#'  and gene expression profiles. This core function of the SigFun package
#'  bridges phenotypic signatures with transcriptomic data to enable downstream
#'  functional analysis.
#'
#' @param seData A `SummarizedExperiment` object containing:
#' - `assays`: Gene expression matrix (genes in rows, samples in columns)
#' - `colData`: Signature scores/classifications for each sample
#' - `rowData`: Gene annotation with ENSEMBL IDs and symbols
#' @param corMethod Correlation method:
#' - "spearman" (default): Non-parametric rank correlation
#' - "pearson": Linear correlation
#' - "kendall": Rank correlation for small datasets
#' - "logit": Logistic regression for binary signatures (e.g., high/low risk)
#' @param zTransform Logical. Whether to z-score normalize expression data:
#' - `FALSE` (default): No transformation (recommended for binary signatures)
#' - `TRUE`: Standardize expression values (recommended for continuous
#' signatures)
#'
#' @return Returns a modified `SummarizedExperiment` object with:
#' - `cor.df` in metadata: Dataframe containing correlation statistics with
#' columns:
#'   - `gene`: ENSEMBL ID
#'   - `cor`: Correlation coefficient (or log-odds for logit)
#'
#' @details This function performs three key operations:
#' 1. Integrates signature data with expression matrices
#' 2. Computes genome-wide correlations using specified method
#' 3. Stores results in SE object for downstream analysis with sig2GSEA
#'
#' @importFrom SummarizedExperiment colData assay rowData
#' @importFrom S4Vectors metadata
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather pivot_wider
#' @importFrom dplyr filter select distinct
#' @importFrom stats cor binomial glm
#' @export
#'
#' @examples
#' data("expr.data")
#' data("mapping")
#' data("SIG_MAT")
#' seData <- SummarizedExperiment::SummarizedExperiment(
#'     assays=list(abundance=as.matrix(expr.data)),
#'     rowData=S4Vectors::DataFrame(mapping, row.names=mapping$ensg_id),
#'     colData=S4Vectors::DataFrame(SIG_MAT))
#'
#' # Binary signature analysis (MammaPrint high/low risk)
#' seCorLogit <- sigCor(seData=seData, corMethod="logit")
#' # Access results
#' S4Vectors::metadata(seCorLogit)$cor.df
sigCor <- function(seData, corMethod="spearman", zTransform=FALSE) {
    signature.obj <- as.data.frame(SummarizedExperiment::colData(seData))

    exp_data <- as.data.frame(SummarizedExperiment::assay(seData))
    exp_data$ensg_id <- rownames(exp_data)
    exp_data <- exp_data %>%
        tidyr::gather(-ensg_id, key="sample_id", value="value")
    S_ID <- intersect(signature.obj$sample_id, unique(exp_data$sample_id))
    exp_data.sid <- exp_data %>% dplyr::filter(sample_id %in% S_ID)
    signature.obj <- signature.obj %>% dplyr::filter(sample_id %in% S_ID)
    exp_data.FPKM_UQ.tmp <- exp_data.sid %>%
        dplyr::select(sample_id, ensg=ensg_id, value)
    tmp.data.wide <- exp_data.FPKM_UQ.tmp %>% dplyr::filter(ensg != "") %>%
        dplyr::select(sample_id, ensg, value) %>%
        dplyr::distinct(sample_id, ensg, .keep_all=TRUE) %>%
        tidyr::pivot_wider(names_from=ensg, values_from=value)

    cor.object <- merge(signature.obj, tmp.data.wide)
    pattern.signature <- cor.object$value
    pattern.genes <- cor.object %>% dplyr::select(-sample_id, -value)
    count.EffectSamples <- unlist(lapply(pattern.genes,
                                         function(x) length(unique(x))))
    rm.index <- which(count.EffectSamples < 2)

    if (zTransform == TRUE) {
        pattern.genes.norm <- if (length(rm.index) > 0) {
            apply(pattern.genes[, -rm.index], 2, .z_score_cal)
        } else { apply(pattern.genes, 2, .z_score_cal) }
    } else { pattern.genes.norm <- if (length(rm.index) > 0)
    {pattern.genes[, -rm.index]} else {pattern.genes} }

    if(corMethod %in% c("pearson", "kendall", "spearman")){
        cor.list <- .corList(pattern.signature, pattern.genes.norm, corMethod)
    }

    if(corMethod %in% c("logit")){
        cor.list <- .logitList(y=pattern.signature, pattern.genes.norm, corMethod)
    }

    cor.df <- data.frame(gene=rownames(cor.list), cor=as.numeric(cor.list$cor) ,
                         pval=as.numeric(cor.list$pval))
    S4Vectors::metadata(seData) <- list(cor.df=cor.df)
    return(seData)
}
