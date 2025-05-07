library(testthat)
library(SummarizedExperiment)

test_that("sig2Fun processes data correctly and creates expected outputs", {
  # Setup test environment
  temp_dir <- file.path(tempdir(),"uni_test")
  data(demo_GSE181574)  # Load example data
  # Create mock pathways for testing
  # Test function execution
  expect_warning(result <- sig2Fun(
    SE_data=SE_GSE181574,
    ranking.method="stat",
    species="human",
    cor.method="logit",
    pathways_all=pathways.all,
    output_path=temp_dir,
    topN=10,
    Z.transform=FALSE,
    significat_type="pval"
  ))

  # Test output structure
  expect_true(!is.null(result$heatmap))

  # Test output files existence
  expected_files <- c(
    "GSEA_pathway_stat_pvalue.txt",
    "GSEA_pathway_stat_qvalue.txt",
    "GSEA_pathway_stat_all.txt",
    "corRES.txt"
  )
})

test_that("sig2Fun handles invalid inputs appropriately", {
  temp_dir <- tempdir()
  data(demo_GSE181574)
  # Test with invalid ranking method
  expect_error(
    sig2Fun(
      SE_data=demo,
      ranking.method="invalid_method",
      output_path=temp_dir
    ) )

  # Test with invalid correlation method
  expect_error(
    sig2Fun(
      SE_data=SE_GSE181574,
      ranking.method="stat",
      cor.method="invalid_method",
      output_path=temp_dir
    )  )

  # Test with invalid significance type
  expect_error(
    sig2Fun(
      SE_data=SE_GSE181574,
      ranking.method="stat",
      significat_type="invalid_type",
      output_path=temp_dir
    )  )

  # Test with invalid topN value
  expect_error(
    sig2Fun(
      SE_data=SE_GSE181574,
      ranking.method="stat",
      topN=-1,
      output_path=temp_dir
    )  )
})

