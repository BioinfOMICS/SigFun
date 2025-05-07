#' @title Example Transcriptomic Dataset for SigFun Demonstration
#' @description A SummarizedExperiment object containing simulated gene expression data
#' and clinical signatures from study GSE180962. This dataset demonstrates SigFun's
#' functionality for analyzing continuous prognostic signatures against transcriptome profiles.
#'
#' @format A `SummarizedExperiment` object with:
#' \describe{
#'   \item{assays}{
#'     - "expression": Gene expression matrix (17,341 genes × 233 samples)
#'       - Rows: Genes (ENSEMBL IDs)
#'       - Columns: Tumor samples
#'       - Values: Normalized log2(CPM) expression values
#'   }
#'   \item{rowData}{
#'     - ensg_id: ENSEMBL gene identifiers
#'     - gene_symbol: Official gene symbols
#'     - gene_biotype: Gene classification (e.g., "protein_coding", "lncRNA")
#'   }
#'   \item{colData}{
#'     - sample_id: Unique sample identifiers
#'     - risk_score: Continuous prognostic risk score (range: 0-1)
#'   }
#'   \item{metadata}{
#'     - Normalization method: RUVseq with 20 factors
#'     - Batch correction: Combat-adjusted
#'   }
#' }
#'
#' @source Synthetic data derived from GSE180962 study patterns,
#' simulating a non-small cell lung cancer cohort. Original study available from GEO:
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180962}
#'
#' @examples
#' # Load dataset
#' data("demo_GSE180962")
#'
#' # Access expression matrix
#' head(SummarizedExperiment::assay(SE_GSE180962)[, 1:3])
#'
#' # View clinical signatures
#' head(SummarizedExperiment::colData(SE_GSE180962))
#'
#' # Check gene annotations
#' head(SummarizedExperiment::rowData(SE_GSE180962))
#'
#' @docType data
#' @keywords datasets
#' @name SE_GSE180962
#' @usage data(demo_GSE180962)
NULL
