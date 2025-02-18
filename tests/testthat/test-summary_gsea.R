library(testthat)
library(dplyr)
library(data.table)

test_that(".summary_gsea function works correctly", {
  # Create a temporary directory for test outputs
  test_output_dir <- tempdir()

  # Create mock GSEA results
  mock_fgsea_res <- data.frame(
    pathway = c("Path1", "Path2", "Path3", "Path4"),
    pval = c(0.01, 0.06, 0.03, 0.002),
    padj = c(0.05, 0.2, 0.1, 0.01),
    NES = c(1.5, -1.2, 0.8, 2.0),
    stringsAsFactors = FALSE
  )

  # Temporarily suppress global assignment warnings
  suppressWarnings({
    # Test with 'stat' ranking method
    expect_no_error({
      .summary_gsea(
        fgseaRes = mock_fgsea_res,
        ranking.method = 'stat',
        output_path = test_output_dir
      )
    })

    # Test with 'abs.stat' ranking method
    expect_no_error({
      .summary_gsea(
        fgseaRes = mock_fgsea_res,
        ranking.method = 'abs.stat',
        output_path = test_output_dir
      )
    })
  })

  # Verify file outputs for 'stat' method
  stat_pvalue_file <- file.path(
    test_output_dir,
    "GSEA_pathway_stat_pvalue.txt"
  )
  stat_qvalue_file <- file.path(
    test_output_dir,
    "GSEA_pathway_stat_qvalue.txt"
  )
  stat_all_file <- file.path(
    test_output_dir,
    "GSEA_pathway_stat_all.txt"
  )

  # Check if files were created
  expect_true(file.exists(stat_pvalue_file))
  expect_true(file.exists(stat_qvalue_file))
  expect_true(file.exists(stat_all_file))

  # Read and verify p-value filtered results
  pvalue_results <- fread(stat_pvalue_file)
  expect_true(nrow(pvalue_results) > 0)
  expect_true(all(pvalue_results$pval < 0.05))

  # Read and verify q-value filtered results
  qvalue_results <- fread(stat_qvalue_file)
  expect_true(nrow(qvalue_results) > 0)
  expect_true(all(qvalue_results$padj < 0.05))

  # Test error handling
  expect_error({
    .summary_gsea(
      fgseaRes = NULL,
      ranking.method = 'invalid',
      output_path = test_output_dir
    )
  })
})
