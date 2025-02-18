library(testthat)
library(SummarizedExperiment)
library(dplyr)
library(tidyverse)
library(fgsea)

# Create dummy data for testing
set.seed(123)  # for reproducibility
data(demo)
# Create a temporary output directory
output_path <- tempdir()

test_that("Sig2GSEA functions correctly", {
  # Setup test data ----
  # 建立模擬的 SummarizedExperiment 物件
  test_matrix <- matrix(rnorm(100), nrow = 10)
  rownames(test_matrix) <- paste0("ENSG", seq_len(10))
  colnames(test_matrix) <- paste0("sample", seq_len(10))

  col_data <- data.frame(
    ensg_id = paste0("ENSG", seq_len(10)),
    gene_symbol = paste0("gene", seq_len(10))
  )

  test_SE <- SummarizedExperiment(
    assays = list(counts = test_matrix),
    colData = col_data
  )

  # 加入模擬的相關性資料
  test_SE@metadata$cor.df <- data.frame(
    gene = paste0("ENSG", seq_len(10)),
    cor = runif(10, -1, 1)
  )

  # 建立測試用的 pathway 資料
  test_pathways <- list(
    pathway1 = paste0("gene", seq_len(3)),
    pathway2 = paste0("gene", 4:6)
  )

  # 建立暫時的輸出目錄
  temp_dir <- tempdir()

  # 執行函數
  result <- Sig2GSEA(
    SE_data.cor = test_SE,
    ranking.method = "stat",
    output_path = temp_dir,
    pathways.all = test_pathways
  )

  # 測試檢查點 ----
  # 1. 檢查回傳值是否為 SummarizedExperiment 物件
  expect_s4_class(result, "SummarizedExperiment")

  # 2. 檢查 metadata 是否包含 fgseaRes
  expect_true("fgseaRes" %in% names(result@metadata))

  # 3. 檢查 fgseaRes 的結構
  expect_true(all(c("pathway", "pval", "padj", "ES", "NES") %in%
                    colnames(result@metadata$fgseaRes)))

  # 4. 檢查是否有產生輸出檔案
  expect_true(file.exists(file.path(temp_dir, "ranks.notdeframe.txt")))

  # 5. 檢查 fgseaRes 的數值是否在合理範圍內
  fgsea_res <- result@metadata$fgseaRes
  expect_true(all(fgsea_res$pval >= 0 & fgsea_res$pval <= 1))
  expect_true(all(fgsea_res$padj >= 0 & fgsea_res$padj <= 1))

  # 6. 測試錯誤處理
  # 測試缺少必要參數時的錯誤
  expect_error(Sig2GSEA())

  # 清理暫時檔案
  unlink(file.path(temp_dir, "ranks.notdeframe.txt"))
})

test_that("Sig2GSEA handles alternative column names", {
  # Setup test data with alternative column names ----
  test_matrix <- matrix(rnorm(100), nrow = 10)
  rownames(test_matrix) <- paste0("ENSG", seq_len(10))
  colnames(test_matrix) <- paste0("sample", seq_len(10))

  col_data <- data.frame(
    V2 = paste0("ENSG", seq_len(10)),
    V1 = paste0("gene", seq_len(10))
  )

  test_SE <- SummarizedExperiment(
    assays = list(counts = test_matrix),
    colData = col_data
  )

  test_SE@metadata$cor.df <- data.frame(
    gene = paste0("ENSG", seq_len(10)),
    cor = runif(10, -1, 1)
  )

  test_pathways <- list(
    pathway1 = paste0("gene", seq_len(3)),
    pathway2 = paste0("gene", 4:6)
  )

  temp_dir <- tempdir()

  # 執行函數
  result <- Sig2GSEA(
    SE_data.cor = test_SE,
    ranking.method = "stat",
    output_path = temp_dir,
    pathways.all = test_pathways
  )

  # 測試檢查點 ----
  expect_s4_class(result, "SummarizedExperiment")
  expect_true("fgseaRes" %in% names(result@metadata))

  # 清理暫時檔案
  unlink(file.path(temp_dir, "ranks.notdeframe.txt"))
})
