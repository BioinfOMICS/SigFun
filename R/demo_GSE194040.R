#' Example SummarizedExperiment Data (GSE194040)
#'
#' @description
#' A `SummarizedExperiment` object containing example gene expression data and sample annotations, provided for demonstrating the functions in the SigFun package. This dataset is structured to facilitate signature-to-function analysis workflows, including correlation, enrichment, and visualization steps.
#'
#' @format A `SummarizedExperiment` object with:
#' \describe{
#'   \item{assays}{A list of matrices containing gene expression data. Each matrix has genes as rows and samples as columns. Use \code{assays(SE_GSE194040)} or \code{assay(SE_GSE194040, "abundance")} to access.}
#'   \item{rowData}{A \code{DataFrame} containing gene-level annotations, typically including:}
#'     \itemize{
#'       \item \code{ensg_id}: ENSEMBL gene identifier (should match rownames of the assay matrix)
#'       \item \code{gene_symbol}: Official gene symbol
#'       \item \code{gene_biotype}: Gene type (e.g., "protein_coding", "lncRNA")
#'     }
#'   \item{colData}{A \code{DataFrame} containing sample-level information, typically including:}
#'     \itemize{
#'       \item \code{sample_id}: Unique sample identifier (should match colnames of the assay matrix)
#'       \item \code{value}: Signature score or classification for each sample (binary or continuous, depending on the analysis)
#'     }
#' }
#'
#' @details
#' This dataset is formatted to be compatible with the SigFun workflow, enabling users to:
#' \itemize{
#'   \item Calculate genome-wide correlations between signatures and gene expression (\code{sigCor})
#'   \item Perform gene set enrichment analysis (\code{sig2GSEA})
#'   \item Visualize results with bar plots and heatmaps (\code{plot_bar}, \code{plot_heat})
#' }
#' The structure mirrors typical transcriptomic datasets used in clinical biomarker research, with coordinated sample and gene annotations to ensure robust downstream analysis[2][1].
#'
#' @source
#' Synthetic data generated for demonstration and testing purposes.
#'
#' @examples
#' # Load the data
#' data(demo_GSE194040)
#'
#' # Explore assay data
#' head(SummarizedExperiment::assay(SE_GSE194040)[, 1:3])
#'
#' # Explore sample annotations
#' head(SummarizedExperiment::colData(SE_GSE194040))
#'
#' # Explore gene annotations
#' head(SummarizedExperiment::rowData(SE_GSE194040))
#'
#' @docType data
#' @keywords datasets
#' @name SE_GSE194040
#' @usage data(demo_GSE194040)
NULL
