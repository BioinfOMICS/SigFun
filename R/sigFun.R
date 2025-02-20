#' SigFun: Functional analysis of Gene Signature
#'
#' SigFun is an R package designed to streamline the analysis of transcriptomic
#' data in relation to specific gene signatures. It provides an automated
#' workflow for analyzing gene signatures and their functional implications in
#' transcriptomic datasets. The package integrates Gene Set Enrichment Analysis
#' (GSEA) with visualization tools to help researchers understand the
#' biological pathways and processes associated with their gene signatures of
#' interest.
#'
#' @section Main functions:
#' - \code{\link{sig2Fun}}: One-click function of running the pipeline
#' Analysis.
#' - \code{\link{sigCor}}: Correlation analysis between signature and genes
#' - \code{\link{sig2GSEA}}: GSEA analysis of signature
#' - \code{\link{plot_bar}}: bar plot visualization of sig2GSEA results
#' - \code{\link{plot_heat}}: heatmap visualization of sig2GSEA results
#'
#' @section Vignettes:
#' See the package vignettes for detailed workflows:
#' \code{vignette('sigFun')}
#'
#' @section Installation:
#' To install from Bioconductor, use:
#' \preformatted{
#' if (!requireNamespace('BiocManager', quietly=TRUE))
#'     install.packages('BiocManager')
#' BiocManager::install('sigFun')
#' }
#'
#' @name sigFun
#' @keywords internal
"_PACKAGE"
