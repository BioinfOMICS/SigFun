library(testthat)
library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(ggplot2)

test_that("plot_heat function produces expected outputs", {

  # Create temporary directory for output
  temp_dir <- tempdir()

  # Create mock pathways
  mock_pathways <- list(
    GOBP_pathway1 = c("gene1", "gene2"),
    GOCC_pathway1 = c("gene3", "gene4"),
    GOMF_pathway1 = c("gene5", "gene6"),
    KEGG_pathway1 = c("gene7", "gene8")
  )
  # Test 1: Basic functionality
  data(demo)
  mock_pathways <- pathways.all

  test_SE <- SE_data
  result <- suppressWarnings(plot_heat(
    SE_data.fgsea = test_SE,
    output_path = temp_dir,
    pathways.all = mock_pathways
  ))

  expect_true(is.list(result))
  expect_equal(length(result), 6)  # Should have 6 elements for default strings

  # Test 2: Check if output files are created
  expected_files <- paste0(
    "GSEA_heatmap_",
    c("KEGG"),
    ".png"
  )

  # Test 3: Test with different significance type
  result_padj <- suppressWarnings(plot_heat(
    SE_data.fgsea = test_SE,
    output_path = temp_dir,
    significat_type = "padj",
    pathways.all = mock_pathways
  ))
  expect_true(is.list(result_padj))

  # Test 4: Test with different topN
  result_top5 <- suppressWarnings(plot_heat(
    SE_data.fgsea = test_SE,
    output_path = temp_dir,
    topN = 5,
    pathways.all = mock_pathways
  ))
  expect_true(is.list(result_top5))

  # Test 5: Test with subset of strings
  result_subset <- plot_heat(
    SE_data.fgsea = test_SE,
    output_path = temp_dir,
    strings = c("KEGG"),
    pathways.all = mock_pathways
  )

  # Test 6: Test error handling for invalid input
  expect_error(
    plot_heat(
      SE_data.fgsea = "invalid_input",
      output_path = temp_dir,
      pathways.all = mock_pathways
    )
  )

  # Clean up temporary files
  unlink(file.path(temp_dir, expected_files))
})
