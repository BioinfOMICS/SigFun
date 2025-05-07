#' @title Breast Cancer Cohort Demo Dataset
#' @description A SummarizedExperiment object containing gene expression data
#' and clinical signatures from the GSE181574 study. This dataset demonstrates
#' SigFun's functionality for analyzing binary clinical signatures (MammaPrint
#' risk classification) against transcriptome data.
#'
#' @format A `SummarizedExperiment` object with:
#' \describe{
#'   \item{assays}{
#'     - "abundance": Gene expression matrix (17,341 genes × 127 samples)
#'       - Rows: Genes (ENSEMBL IDs)
#'       - Columns: Tumor samples
#'       - Values: Normalized expression values
#'   }
#'   \item{rowData}{
#'     - ensg_id: ENSEMBL gene identifiers
#'     - gene_symbol: Official gene symbols
#'     - gene_biotype: Gene classification (e.g., "protein_coding", "lncRNA")
#'   }
#'   \item{colData}{
#'     - sample_id: Unique sample identifiers
#'     - value: Binary risk classification (0=low risk, 1=high risk)
#'   }
#'   \item{metadata}{
#'     - Study design information
#'     - Signature interpretation parameters
#'   }
#' }
#'
#' @source Synthetic data derived from GSE181574 study (PMID: 34903889),
#' modified for demonstration purposes. Original data available from GEO:
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181574}
#'
#' @examples
#' # Load dataset
#' data("demo_GSE181574")
#'
#' # Access expression matrix
#' head(SummarizedExperiment::assay(SE_GSE181574)[, 1:3])
#'
#' # View clinical signatures
#' head(SummarizedExperiment::colData(SE_GSE181574))
#'
#' # Check gene annotations
#' head(SummarizedExperiment::rowData(SE_GSE181574))
#'
#' @docType data
#' @keywords datasets
#' @name SE_GSE181574
#' @usage data(demo_GSE181574)
NULL
