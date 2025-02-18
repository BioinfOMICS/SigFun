library(testthat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)

# Create sample test data
create_test_data <- function() {
  data.frame(
    pathway = paste0("pathway_", seq_len(15)),
    NES = c(rep(1.5, 8), rep(-1.5, 7)),
    pvalue = runif(15),
    stringsAsFactors = FALSE
  )
}

# Test suite for .barplot function
test_that(".barplot processes data and creates plot correctly", {
  # Setup
  temp_dir <- tempdir()
  test_data <- create_test_data()
  plot_name <- "test_plot.pdf"

  # Test execution
  .barplot(
    plot.name = plot_name,
    type.sig = "pvalue",
    topN = 5,
    RES_NES_strings = test_data,
    output_path = temp_dir
  )

  # Test file creation
  expect_true(
    file.exists(file.path(temp_dir, plot_name)),
    "Plot file should be created"
  )

  # Test data processing
  processed_data <- test_data %>%
    dplyr::filter(NES > 0) %>%
    dplyr::slice(seq_len(5))

  #expect_equal(
  #  nrow(processed_data),
  #  5,
  #  "Should return top 5 positive NES values"
  #)

  # Clean up
  unlink(file.path(temp_dir, plot_name))
})

test_that(".barplot handles edge cases correctly", {
  # Setup
  temp_dir <- tempdir()

  # Test with empty data
  empty_data <- data.frame(
    pathway = character(),
    NES = numeric(),
    pvalue = numeric(),
    stringsAsFactors = FALSE
  )

  # Should not error with empty data
  #expect_error(
  #  .barplot(
  #    plot.name = "empty_plot.pdf",
  #    type.sig = "pvalue",
  #    topN = 5,
  #    RES_NES_strings = empty_data,
  #    output_path = temp_dir
  #  ),
  #  NA
  #)

  # Test with single row
  single_row_data <- data.frame(
    pathway = "single_pathway",
    NES = 1.5,
    pvalue = 0.05,
    stringsAsFactors = FALSE
  )

  expect_error(
    .barplot(
      plot.name = "single_plot.pdf",
      type.sig = "pvalue",
      topN = 5,
      RES_NES_strings = single_row_data,
      output_path = temp_dir
    ),
    NA
  )
})

test_that(".barplot respects topN parameter", {
  # Setup
  temp_dir <- tempdir()
  test_data <- create_test_data()

  # Test with different topN values
  test_cases <- c(3, 5, 10, 20)

  for (n in test_cases) {
    plot_name <- paste0("test_plot_", n, ".pdf")

    .barplot(
      plot.name = plot_name,
      type.sig = "pvalue",
      topN = n,
      RES_NES_strings = test_data,
      output_path = temp_dir
    )

    # Verify file exists
    expect_true(
      file.exists(file.path(temp_dir, plot_name)),
      paste("Plot file should be created for topN =", n)
    )

    # Clean up
    unlink(file.path(temp_dir, plot_name))
  }
})
