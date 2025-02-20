library(testthat)
library(dplyr)
library(ggplot2)
library(cowplot)

test_that(".plotGseaTable handles inputs correctly and returns expected output", {
  # Setup test data
  set.seed(42)
  stats <- setNames(rnorm(100), paste0("gene", seq_len(100)))
  pathways <- list(
    "pathway1"=paste0("gene", seq_len(10)),
    "pathway2"=paste0("gene", seq(20, 30))  # Fixed: using seq() instead of seq_len()
  )

  fgseaRes <- data.frame(
    pathway=names(pathways),
    pval=c(0.001, 0.002),
    padj=c(0.01, 0.02),
    NES=c(1.5, -2.0),
    size=c(10, 11)
  )

  # Test 1: Basic functionality
  result <- .plotGseaTable(
    pathways_all=pathways,
    Pathways=pathways,
    stats=stats,
    fgseaRes=fgseaRes
  )

  expect_true(inherits(result, "ggplot"))

  # Test 2: Custom style parameters
  custom_pathway_style <- list(size=14, hjust=0.1, x=0.1, vjust=0.6)
  result_custom <- .plotGseaTable(
    pathways_all=pathways,
    Pathways=pathways,
    stats=stats,
    fgseaRes=fgseaRes,
    pathwayLabelStyle=custom_pathway_style
  )

  expect_true(inherits(result_custom, "ggplot"))


  # Test 5: Custom column widths
  custom_widths <- c(1, 1, 1.2, 3.5, 12)
  result_widths <- .plotGseaTable(
    pathways_all=pathways,
    Pathways=pathways,
    stats=stats,
    fgseaRes=fgseaRes,
    colwidths=custom_widths
  )

  expect_true(inherits(result_widths, "ggplot"))

  # Test 6: Deprecated render parameter
  expect_warning(
    .plotGseaTable(
      pathways_all=pathways,
      Pathways=pathways,
      stats=stats,
      fgseaRes=fgseaRes,
      render=TRUE
    ),
    "render argument is deprecated"
  )

  # Test 7: gseaParam effect
  result_param <- .plotGseaTable(
    pathways_all=pathways,
    Pathways=pathways,
    stats=stats,
    fgseaRes=fgseaRes,
    gseaParam=2
  )

  expect_true(inherits(result_param, "ggplot"))

})

# Helper function to create mock data for testing
create_mock_gsea_data <- function(n_genes=100, n_pathways=2) {
  stats <- setNames(rnorm(n_genes), paste0("gene", seq_len(n_genes)))

  pathways <- lapply(seq_len(n_pathways), function(i) {
    paste0("gene", seq_len(sample(10:20, 1)))
  })
  names(pathways) <- paste0("pathway", seq_len(n_pathways))

  fgseaRes <- data.frame(
    pathway=names(pathways),
    pval=runif(n_pathways, 0, 0.05),
    padj=runif(n_pathways, 0, 0.1),
    NES=rnorm(n_pathways, 0, 2),
    size=vapply(pathways, length, numeric(1))
  )

  list(
    stats=stats,
    pathways=pathways,
    fgseaRes=fgseaRes
  )
}
