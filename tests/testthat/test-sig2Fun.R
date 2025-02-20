library(testthat)
library(SummarizedExperiment)

test_that("sig2Fun processes data correctly and creates expected outputs", {
  # Setup test environment
  temp_dir <- file.path(tempdir(),"uni_test")
  data(demo)  # Load example data
  # Create mock pathways for testing
  SE_data@metadata$cor.df <- NULL
  SE_data@metadata$fgseaRes <- NULL
  # Test function execution
  expect_warning(result <- sig2Fun(
    SE_data=SE_data,
    ranking.method="stat",
    species="human",
    cor.method="spearman",
    pathways_all=pathways.all,
    output_path=temp_dir,
    topN=10,
    Z.transform=FALSE,
    significat_type="pval",
    strings=c("KEGG"),
    plot_out=TRUE
  ))

  # Test output structure
  expect_true(!is.null(result$heatmap))
  expect_true(!is.null(result@metadata$fgseaRes))
  expect_true(!is.null(result@metadata$cor.df))

  # Test output directory creation and permissions
  expect_true(dir.exists(temp_dir))

  # Test output files existence
  expected_files <- c(
    "GSEA_pathway_stat_pvalue.txt",
    "GSEA_pathway_stat_qvalue.txt",
    "GSEA_pathway_stat_all.txt",
    "corRES.txt"
  )

  for (file in expected_files) {
    expect_true(
      file.exists(file.path(temp_dir, file)),
      info=paste("Missing expected output file:", file)
    )
  }
})

test_that("sig2Fun handles invalid inputs appropriately", {
  temp_dir <- tempdir()
  data(demo)
  SE_data@metadata$cor.df <- NULL
  SE_data@metadata$fgseaRes <- NULL
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
      SE_data=demo,
      ranking.method="stat",
      cor.method="invalid_method",
      output_path=temp_dir
    )  )

  # Test with invalid significance type
  expect_error(
    sig2Fun(
      SE_data=demo,
      ranking.method="stat",
      significat_type="invalid_type",
      output_path=temp_dir
    )  )

  # Test with invalid topN value
  expect_error(
    sig2Fun(
      SE_data=demo,
      ranking.method="stat",
      topN=-1,
      output_path=temp_dir
    )  )
})

test_that("sig2Fun handles parameter variations correctly", {
  temp_dir <- tempdir()
  data(demo)
  SE_data@metadata$cor.df <- NULL
  SE_data@metadata$fgseaRes <- NULL
  mock_pathways <- list(
    GOBP_PATH1=c("GENE1", "GENE2", "GENE3")
  )

  # Test with different correlation methods
  expect_warning(result_pearson <- sig2Fun(
    SE_data=SE_data,
    ranking.method="stat",
    cor.method="pearson",
    pathways_all=mock_pathways,
    output_path=temp_dir,
    Z.transform=TRUE
  ))

  expect_warning(result_spearman <- sig2Fun(
    SE_data=SE_data,
    ranking.method="stat",
    cor.method="spearman",
    pathways_all=mock_pathways,
    output_path=temp_dir,
    Z.transform=TRUE
  ))

  expect_true(
      sum(abs(result_pearson@metadata$cor.df$cor-
          result_spearman@metadata$cor.df$cor)) != 0
  )

  # Clean up
  unlink(temp_dir, recursive=TRUE)
})
