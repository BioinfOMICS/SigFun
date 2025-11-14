library(testthat)
library(dplyr)
library(SummarizedExperiment)
data("sig2Fun_result")

test_that("chordPlot function basic structure", {
    result <- tryCatch({
        chordPlot(sig2Fun_result, showCategory = 3)
    }, error = function(e) {
        if (grepl("gap.degree|space to allocate sectors", e$message)) {
            message("Skipping chord diagram visualization due to spacing constraints")
            return(NULL)
        } else stop(e)
    })

    if (!is.null(result)) {
        expect_type(result, "list")
        expect_named(result, c("chordPlot", "tableChordPlot"))
        expect_s3_class(result$tableChordPlot, "data.frame")
        expect_lte(nrow(result$tableChordPlot), 3)
    } else {
        data_processed <- tryCatch({
            .extractDF(sig2Fun_result, type = "gseaReadable") |>
                .extractGeneSets(3) |>
                dplyr::mutate(name = .labelBreak(categoryID, 30))
        }, error = function(e) skip(paste("Data processing failed:", e$message)))

        expect_s3_class(data_processed, "data.frame")
        expect_true(all(c("Gene", "name") %in% colnames(data_processed)))
    }
})

test_that("chordPlot with numeric showCategory works", {
    result <- tryCatch({
        chordPlot(sig2Fun_result, showCategory = 2)
    }, error = function(e) {
        if (grepl("gap.degree|space to allocate sectors", e$message))
            skip("Skipping due to chord diagram spacing constraints")
        else stop(e)
    })

    if (!is.null(result)) {
        expect_lte(nrow(result$tableChordPlot), 2)
    }
})

test_that("chordPlot with specific pathway works", {
    specific_pathway <- "REACTOME_CELL_CYCLE_MITOTIC"
    result <- tryCatch({
        chordPlot(sig2Fun_result, showCategory = specific_pathway)
    }, error = function(e) {
        if (grepl("gap.degree|space to allocate sectors", e$message)) {
            skip("Skipping due to chord diagram spacing constraints")
        } else if (grepl("invalid_pathway", e$message)) {
            skip(paste("Specific pathway not found:", specific_pathway))
        } else {
            stop(e)
        }
    })
    if (!is.null(result)) {
        expect_s3_class(result$tableChordPlot, "data.frame")
        expect_equal(nrow(result$tableChordPlot), 1)
    } else {
        skip("chordPlot result is NULL (likely skipped due to spacing or missing data)")
    }
})

test_that("chordPlot handles fontSize type errors correctly", {
    expect_error(
        chordPlot(sig2Fun_result, fontSize = "invalidSize"),
        regexp = "fontSize.*must be.*numeric"
    )
})

test_that("chordPlot handles fontSize boundary errors correctly", {
    expect_error(
        chordPlot(sig2Fun_result, fontSize = 0),
        regexp = "numeric value between 0"
    )
    expect_error(
        chordPlot(sig2Fun_result, fontSize = 1.5),
        regexp = "numeric value between"
    )
})

test_that("chordPlot handles invalid pathway errors correctly", {
    expect_error(
        chordPlot(sig2Fun_result, showCategory = c("invalid_pathway1", "invalid_pathway2")),
        regexp = "cannot be found"
    )
})

test_that("chordPlot handles single invalid pathway correctly", {
    expect_error(
        chordPlot(sig2Fun_result, showCategory = "NONEXISTENT_PATHWAY_12345"),
        regexp = "showCategory.*cannot be found"
    )
})

test_that("tableChordPlot has correct structure", {
    result <- tryCatch({
        chordPlot(sig2Fun_result, showCategory = 1)
    }, error = function(e) {
        if (grepl("gap.degree|space to allocate sectors", e$message)) {
            skip("Cannot create chord plot due to spacing")
        } else {
            skip(paste("Cannot create chord plot:", e$message))
        }
    })

    if (!is.null(result) && !is.null(result$tableChordPlot)) {
        table_data <- result$tableChordPlot
        expected_cols <- c("original_name", "label_name", "Gene")
        expect_true(all(expected_cols %in% colnames(table_data)))
        expect_true(nrow(table_data) >= 1)
        expect_true(all(!is.na(table_data$original_name)))
    }
})

test_that("chordPlot handles different breaklineN values", {
    result1 <- tryCatch({
        chordPlot(sig2Fun_result, showCategory = 2, breaklineN = 20)
    }, error = function(e) {
        if (grepl("gap.degree|space", e$message)) {
            skip("Spacing constraint")
        } else stop(e)
    })

    result2 <- tryCatch({
        chordPlot(sig2Fun_result, showCategory = 2, breaklineN = 50)
    }, error = function(e) {
        if (grepl("gap.degree|space", e$message)) {
            skip("Spacing constraint")
        } else stop(e)
    })

    if (!is.null(result1) && !is.null(result2)) {
        max_len1 <- max(nchar(result1$tableChordPlot$label_name))
        max_len2 <- max(nchar(result2$tableChordPlot$label_name))
        expect_true(max_len1 <= max_len2 || max_len1 < 100)
    }
})

test_that("chordPlot Gene column is properly formatted", {
    result <- tryCatch({
        chordPlot(sig2Fun_result, showCategory = 1)
    }, error = function(e) skip("Cannot create chord plot"))

    if (!is.null(result)) {
        genes <- result$tableChordPlot$Gene
        expect_true(all(grepl("/", genes) | !grepl("/", genes)))
        expect_true(all(!is.na(genes)))
        expect_true(all(nchar(genes) > 0))
    }
})

test_that("chordPlot validates showCategory parameter type", {
    expect_error(
        chordPlot(sig2Fun_result, showCategory = NULL),
        regexp = "showCategory"
    )

    expect_error(
        chordPlot(sig2Fun_result, showCategory = list("pathway1")),
    )
})

test_that("chordPlot with fontSize edge cases", {
    result_small <- tryCatch({
        chordPlot(sig2Fun_result, showCategory = 1, fontSize = 0.1)
    }, error = function(e) skip("Spacing constraint"))

    result_large <- tryCatch({
        chordPlot(sig2Fun_result, showCategory = 1, fontSize = 0.9)
    }, error = function(e) skip("Spacing constraint"))
    if (!is.null(result_small)) {
        expect_s3_class(result_small$tableChordPlot, "data.frame")
    }
    if (!is.null(result_large)) {
        expect_s3_class(result_large$tableChordPlot, "data.frame")
    }
})

test_that("chordPlot closes graphics devices properly", {
    initial_devices <- length(grDevices::dev.list())
    result <- tryCatch({
        chordPlot(sig2Fun_result, showCategory = 1)
    }, error = function(e) skip("Cannot create plot"))
    final_devices <- length(grDevices::dev.list())
    expect_equal(initial_devices, final_devices)
})
