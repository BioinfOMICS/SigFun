#' Expression Data
#'
#' @name expr.data
#' @aliases expr.data
#' @docType data
#' @description Expression data for analysis
#' @usage data("expr.data")
#' @format A data frame with 17341 rows and 127 columns
#' @examples data(expr.data)
#' @keywords datasets
NULL

#' Signature Matrix
#'
#' @name SIG_MAT
#' @aliases SIG_MAT
#' @docType data
#' @description A toy dataset representing signature values for samples.
#' @usage data("SIG_MAT")
#' @format A data frame with 127 rows and 2 columns
#' @examples data(SIG_MAT)
#' @keywords datasets
NULL

#' Gene Mapping Table
#'
#' @name mapping
#' @aliases mapping
#' @docType data
#' @description Mapping table for gene identifiers and annotations.
#' @usage data("mapping")
#' @format A data frame with 17341 rows and 3 columns
#' \describe{
#'   \item{ensg_id}{Ensembl gene ID, e.g., ENSG00000123456}
#'   \item{gene_symbol}{Gene symbol, e.g., TP53}
#'   \item{gene_biotype}{Gene biotype, e.g., protein_coding, lncRNA}
#' }
#' @examples data(mapping)
#' @keywords datasets
NULL

#' Transcript to Gene Mapping
#'
#' @name t2g
#' @aliases t2g
#' @docType data
#' @description Another mapping object for transcript to gene relationships.
#' @usage data("t2g")
#' @format A data frame with 17341 rows and 3 columns
#' \describe{
#'   \item{gs_name}{Gene set name or pathway name}
#'   \item{ensembl_gene}{Ensembl gene ID, e.g., ENSG00000123456}
#' }
#' @examples data(t2g)
#' @keywords datasets
NULL

#' SigFun result
#'
#' @name sig2Fun_result
#' @aliases sig2Fun_result
#' @docType data
#' @description A SummarizedExperiment object containing gene expression data
#' and clinical signatures from the GSE181574 study. This dataset demonstrates
#' SigFun's functionality for analyzing binary clinical signatures (MammaPrint
#' risk classification) against transcriptome data.
#' @usage data("sig2Fun_result")
#' @format A SummarizedExperiment  with assay, colData, and metadata.
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
#' @examples data(sig2Fun_result)
#' @keywords datasets
NULL

