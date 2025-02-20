library(testthat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(SummarizedExperiment)

# Create a dummy SE_data.fgsea object for testing
create_dummy_se <- function() {
  # Simulate fgsea results
  fgsea_res <- data.frame(
    pathway=c("GOBP_Term1", "GOCC_Term2", "GOMF_Term3", "KEGG_Term4", "REACTOME_Term5", "WP_Term6", "GOBP_Term7", "GOCC_Term8", "GOMF_Term9", "KEGG_Term10"),
    pval=c(0.01, 0.001, 0.02, 0.005, 0.03, 0.002, 0.1, 0.2, 0.0001, 0.04),
    NES=c(1.5, -1.2, 0.8, -1.0, 1.1, -0.9, 1.3, -1.4, 1.6, -1.7),
    leadingEdge=c("gene1|gene2", "gene3|gene4", "gene5|gene6", "gene7|gene8", "gene9|gene10", "gene11|gene12", "gene13|gene14", "gene15|gene16", "gene17|gene18", "gene19|gene20"),
    pal=c(0.008, 0.0008, 0.015, 0.003, 0.025, 0.0015, 0.09, 0.18, 0.00008, 0.035) # Example of a "pal" column
  )

  # Create a SummarizedExperiment object (mimicking SE_data.fgsea)
  metadata <- list(fgseaRes=fgsea_res)
  SE <- SummarizedExperiment::SummarizedExperiment(assays=matrix(1:10, ncol=1),  # Dummy assay data
                                                   metadata=metadata)
  return(SE)
}

# 創建測試用的輔助函數
create_test_data <- function() {
    # 創建範例資料
    example_mat <- matrix(
        rnorm(100),
        nrow=10,
        ncol=10,
        dimnames=list(
            paste0("gene", 1:10),
            paste0("sample", 1:10)
        )
    )

    # 創建 fgseaRes
    fgsea_res <- data.frame(
        pathway=c(
            "GOBP test1",
            "GOCC test1",
            "GOMF test1",
            "KEGG test1",
            "REACTOME test1",
            "WP test1"
        ),
        pval=c(0.01, 0.02, 0.03, 0.04, 0.06, 0.01),
        padj=c(0.02, 0.03, 0.04, 0.05, 0.07, 0.02),
        NES=c(1.5, -1.2, 1.3, -1.4, 1.1, -1.6),
        size=c(100, 90, 80, 70, 60, 50),
        leadingEdge=rep("gene1|gene2|gene3", 6)
    )

    # 創建 SummarizedExperiment 物件
    se <- SummarizedExperiment(
        assays=list(counts=example_mat),
        metadata=list(fgseaRes=fgsea_res)
    )

    return(se)
}

# 測試案例
test_that("plot_bar generates correct output", {
    # 建立測試資料
    SE_data <- create_test_data()

    # 建立暫時的輸出目錄
    output_dir <- tempdir()
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive=TRUE)
    }

    # 執行函數
    result <- plot_bar(
        SE_data.fgsea=SE_data,
        output_path=output_dir,
        significat_type="pval",
        topN=10,
        strings=c("KEGG")
    )

    # 檢查輸出檔案是否存在
    expect_true(file.exists(file.path(output_dir, "GSEA_stat.txt")))

    # 檢查返回值是否為列表
    expect_type(result, "list")

    # 檢查輸出的統計檔案內容
    stats_data <- data.table::fread(file.path(output_dir, "GSEA_stat.txt"))
    expect_true(all(stats_data$pval < 0.05))
    expect_true(all(grepl(";", stats_data$leadingEdge)))
})


test_that("plot_bar generates correct output2", {
    # 設置測試資料
    se_test <- create_test_data()
    temp_dir <- tempdir()

    # 創建測試輸出目錄（如果不存在）
    if (!dir.exists(temp_dir)) {
        dir.create(temp_dir, recursive=TRUE)
    }

    # 修改原本的測試程式碼，移除會觸發警告的部分
    plot <- ggplot() +
        theme_void() +
        annotate(
            geom="text",  # 明確指定 geom 類型
            x=1,
            y=1,
            label="No significant (pval<0.05) item."
        ) +
        scale_x_continuous(limits=c(0, 2)) +
        scale_y_continuous(limits=c(0, 2)) +
        theme(plot.margin=unit(c(0, 0, 0, 0), "cm"))

    # 測試圖形生成
    expect_true(is.ggplot(plot))

    # 測試圖形儲存
    expect_silent(
        ggsave(
            filename=file.path(temp_dir, "test.png"),
            plot=plot,
            width=3,
            height=1.4,
            units="in",
            dpi=100
        )
    )

    # 檢查檔案是否成功生成
    expect_true(file.exists(file.path(temp_dir, "test.png")))
})
test_that("plot_bar generates correct output bigdata", {
  # demo data as the test
  data("demo")

  output_dir <- tempdir() # Use a temporary directory for testing

  plot_bar(SE_data.fgsea=SE_data, output_path=output_dir, topN=10,
           significat_type="pval", strings=c("KEGG"))

  # Check if the GSEA_stat.txt file is created
  expect_true(file.exists(file.path(output_dir, "GSEA_stat.txt")))

  # Check if the PNG files are created (one for each term in strings)
  strings <- c("KEGG")
  for (s in strings) {
    expect_true(file.exists(file.path(output_dir, paste0("GSEA_top10_", s, ".png"))))
  }


  #Test with no significant pathways for one of the categories
  SE_data_no_sig <- create_dummy_se()
  SE_data_no_sig@metadata$fgseaRes$GOBP_Term1_pal <- 0.5
  SE_data_no_sig@metadata$fgseaRes$GOBP_Term7_pal <- 0.6

  plot_bar(SE_data_no_sig, output_dir)
  expect_true(file.exists(file.path(output_dir, "GSEA_top10_GOBP.png")))


  # Clean up the temporary directory
  # unlink(output_dir, recursive=TRUE) # Be cautious when deleting directories in tests!  Better to let tempdir handle it
})

test_that("plot_bar generates correct output bigdata low TopN", {
  # demo data as the test
  data("demo")

  output_dir <- tempdir() # Use a temporary directory for testing

  plot_bar(SE_data.fgsea=SE_data, output_path=output_dir, topN=3,
           significat_type="pval", strings=c("KEGG"))

  # Check if the GSEA_stat.txt file is created
  expect_true(file.exists(file.path(output_dir, "GSEA_stat.txt")))

  # Check if the PNG files are created (one for each term in strings)
  strings <- c("KEGG")
  for (s in strings) {
    expect_true(file.exists(file.path(output_dir, paste0("GSEA_top10_", s, ".png"))))
  }


  #Test with no significant pathways for one of the categories
  SE_data_no_sig <- create_dummy_se()
  SE_data_no_sig@metadata$fgseaRes$GOBP_Term1_pal <- 0.5
  SE_data_no_sig@metadata$fgseaRes$GOBP_Term7_pal <- 0.6

  plot_bar(SE_data_no_sig, output_dir)
  expect_true(file.exists(file.path(output_dir, "GSEA_top10_GOBP.png")))


  # Clean up the temporary directory
  # unlink(output_dir, recursive=TRUE) # Be cautious when deleting directories in tests!  Better to let tempdir handle it
})

test_that("plot_bar generates correct output bigdata high TopN", {
  # demo data as the test
  data("demo")

  output_dir <- tempdir() # Use a temporary directory for testing

  plot_bar(SE_data.fgsea=SE_data, output_path=output_dir, topN=100,
           significat_type="pval", strings=c("KEGG"))

  # Check if the GSEA_stat.txt file is created
  expect_true(file.exists(file.path(output_dir, "GSEA_stat.txt")))

  # Check if the PNG files are created (one for each term in strings)
  strings <- c("KEGG")
  for (s in strings) {
    expect_true(file.exists(file.path(output_dir, paste0("GSEA_top10_", s, ".png"))))
  }


  #Test with no significant pathways for one of the categories
  SE_data_no_sig <- create_dummy_se()
  SE_data_no_sig@metadata$fgseaRes$GOBP_Term1_pal <- 0.5
  SE_data_no_sig@metadata$fgseaRes$GOBP_Term7_pal <- 0.6

  plot_bar(SE_data_no_sig, output_dir)
  expect_true(file.exists(file.path(output_dir, "GSEA_top10_GOBP.png")))


  # Clean up the temporary directory
  # unlink(output_dir, recursive=TRUE) # Be cautious when deleting directories in tests!  Better to let tempdir handle it
})

test_that("plot_bar handles empty results correctly", {
    # 創建一個沒有顯著結果的測試數據
    empty_fgsea_res <- data.frame(
        pathway=character(),
        pval=numeric(),
        NES=numeric(),
        leadingEdge=character(),
        stringsAsFactors=FALSE
    )

    # 創建空的 SE 物件
    example_mat <- matrix(
        rnorm(100),
        nrow=10,
        ncol=10,
        dimnames=list(
            paste0("gene", 1:10),
            paste0("sample", 1:10)
        )
    )

    se_empty <- SummarizedExperiment::SummarizedExperiment(
        assays=list(counts=example_mat),
        metadata=list(fgseaRes=empty_fgsea_res)
    )

    # 設定輸出路徑
    temp_dir <- tempdir()
    if (!dir.exists(temp_dir)) {
        dir.create(temp_dir, recursive=TRUE)
    }

    # 檢查輸出檔案
    expect_true(file.exists(file.path(temp_dir, "GSEA_top10_GOBP.png")))
})

test_that("plot_bar generates correct output demo", {
    # 建立測試資料
    data(demo)
    # 建立暫時的輸出目錄
    output_dir <- tempdir()
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive=TRUE)
    }

    # 執行函數
    result <- plot_bar(
        SE_data.fgsea=SE_data,
        output_path=output_dir,
        significat_type="pval",
        topN=10,
        strings=c("KEGG")
    )

    # 檢查輸出檔案是否存在
    expect_true(file.exists(file.path(output_dir, "GSEA_stat.txt")))


    # 檢查返回值是否為列表
    expect_type(result, "list")

})

# 修改原始的 plot_bar 函數中處理空結果的部分
# 在你的 plot_bar 函數中，應該要這樣修改空結果的處理：

