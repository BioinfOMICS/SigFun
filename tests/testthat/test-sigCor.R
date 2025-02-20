library(testthat)
library(SummarizedExperiment)
library(tibble)
library(dplyr)
library(tidyr)
library(readr)

test_that("sigCor function works correctly", {
  # Create test data
  # 1. Create a small test SummarizedExperiment object
  test_counts <- matrix(
    rnorm(20),
    nrow=2,  # 2 samples
    ncol=10, # 10 genes
    dimnames=list(
      c("sample1", "sample2"),
      paste0("gene", seq_len(10))
    )
  )



  test_counts3 <- rbind(test_counts, test_counts)
  row.names(test_counts3) <- c("sample1", "sample2", "sample13", "sample4")

  test_rowData <- data.frame(
    sample_id=c("sample1", "sample2"),
    value=c(1.5, 2.5),
    row.names=c("sample1", "sample2")
  )


  test_rowData3 <- data.frame(
    sample_id=c("sample1", "sample2", "sample13", "sample4"),
    value=c(1.5, 2.5, 1.5, 2.5),
    row.names=c("sample1", "sample2", "sample13", "sample4")
  )

  test_SE <- SummarizedExperiment(
    assays=list(counts=test_counts),
    rowData=test_rowData
  )


  test_SE3 <- SummarizedExperiment(
    assays=list(counts=test_counts3),
    rowData=test_rowData3
  )

  # 修正: 直接使用臨時目錄作為輸出路徑
  temp_dir <- tempdir()

  # Test 1: Basic functionality
  result <- sigCor(
    SE_data=test_SE,
    output_path=temp_dir,  # 只傳入目錄路徑
    cor.method="spearman",
    Z.transform=FALSE
  )

  result <- sigCor(
    SE_data=test_SE,
    output_path=temp_dir,  # 只傳入目錄路徑
    cor.method="spearman",
    Z.transform=TRUE
  )


  result <- sigCor(
    SE_data=test_SE3,
    output_path=temp_dir,  # 只傳入目錄路徑
    cor.method="spearman",
    Z.transform=FALSE
  )

  result <- sigCor(
    SE_data=test_SE3,
    output_path=temp_dir,  # 只傳入目錄路徑
    cor.method="spearman",
    Z.transform=TRUE
  )

  data(demo)
  result <- sigCor(
    SE_data=SE_data,
    output_path=temp_dir,  # 只傳入目錄路徑
    cor.method="spearman",
    Z.transform=FALSE
  )

  result <- sigCor(
    SE_data=SE_data,
    output_path=temp_dir,  # 只傳入目錄路徑
    cor.method="spearman",
    Z.transform=TRUE
  )

  # 修正: 更新檢查的檔案路徑
  output_file <- file.path(temp_dir, "corRES.txt")

  # Check if function returns a SummarizedExperiment object
  expect_s4_class(result, "SummarizedExperiment")

  # Check if correlation results exist in metadata
  expect_true("cor.df" %in% names(result@metadata))

  # Check if output file was created
  expect_true(file.exists(output_file))

  # Test 2: Check correlation dataframe structure
  cor_df <- result@metadata$cor.df
  expect_true(all(c("gene", "cor") %in% colnames(cor_df)))
  expect_true(is.numeric(cor_df$cor))

  # Test 3: Test with Z-transform
  result_z <- sigCor(
    SE_data=test_SE,
    output_path=temp_dir,  # 只傳入目錄路徑
    cor.method="spearman",
    Z.transform=TRUE
  )
  expect_s4_class(result_z, "SummarizedExperiment")

  # Test 4: Check error handling for invalid input
  expect_error(
    sigCor(
      SE_data="invalid_input",
      output_path=temp_dir
    )
  )

  # Test 5: Check correlation values are between -1 and 1
  cor_values <- result@metadata$cor.df$cor
  expect_true(all(cor_values >= -1 & cor_values <= 1))

  # Clean up
  unlink(output_file)
})
