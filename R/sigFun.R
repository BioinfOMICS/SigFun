#' SigFun: A Functional Analysis Tool for Signature.
#'
#' SigFun provides a systematic workflow for analyzing biological functions associated with multi-gene signatures, overcoming limitations of conventional methods for small gene sets. The package enables functional interpretation of signatures in various formats (binary classifications, continuous scores) by leveraging whole transcriptome data as a surrogate.
#'
#' @section Key Features:
#' - Analyzes signatures of any size (from single genes to large panels)
#' - Handles both binary (e.g., high/low risk) and continuous signatures
#' - Integrates correlation analysis with pathway enrichment (GSEA)
#' - Generates publication-quality visualizations
#' - Compatible with standard Bioconductor objects (SummarizedExperiment)
#'
#' @section Main Workflow:
#' \enumerate{
#'   \item \code{\link{sigCor}}: Calculate genome-wide correlations between signature and transcriptome
#'   \item \code{\link{sig2GSEA}}: Perform pathway enrichment analysis using correlation statistics
#'   \item \code{\link{plot_bar}}/\code{\link{plot_heat}}: Visualize significant pathways
#'   \item \code{\link{sig2Fun}}: Complete analysis pipeline in one function
#' }
#'
#' @section Vignettes:
#' Detailed workflow and case studies available in:
#' \code{browseVignettes(package = "SigFun")}
#'
#' @section Installation:
#' Install from Bioconductor:
#' \preformatted{
#' if (!requireNamespace("BiocManager", quietly = TRUE))
#'     install.packages("BiocManager")
#' BiocManager::install("SigFun")
#' }
#'
#' @section Datasets:
#' \describe{
#'   \item{\code{SE_GSE181574}}{Breast cancer cohort with MammaPrint classification}
#'   \item{\code{pathways.all}}{MSigDB pathway collection (subset)}
#' }
#'
#' @references
#' For methodology details see the package vignette and:
#' Subramanian A, et al. (2005) Gene set enrichment analysis. PNAS 102:15545-50
#'
#' @name SigFun
#' @docType _PACKAGE
#' @keywords package
NULL
