library(testthat)

testthat::test_that(".defaultLabeller works as expected", {
    labeller <- .defaultLabeller(n = 10)

    testthat::expect_s3_class(labeller, "labeller")
    testthat::expect_type(labeller, "closure")

    result <- labeller("hello_world")
    testthat::expect_true(grepl("hello", result))
    testthat::expect_true(grepl("world", result))
    testthat::expect_false(grepl("_", result, fixed = TRUE))

    result <- labeller("one_two_three")
    testthat::expect_false(grepl("_", result, fixed = TRUE))
    testthat::expect_true(grepl("one", result))
    testthat::expect_true(grepl("two", result))
    testthat::expect_true(grepl("three", result))

    result <- labeller("hello")
    testthat::expect_true(grepl("hello", result))
})

testthat::test_that(".defaultLabeller wraps long strings", {
    labeller <- .defaultLabeller(n = 10)
    long_str <- "very_long_string_that_exceeds_limit"
    result <- labeller(long_str)

    testthat::expect_type(result, "character")
    testthat::expect_false(grepl("_", result, fixed = TRUE))
})

testthat::test_that(".defaultLabeller handles edge cases", {
    labeller <- .defaultLabeller(n = 5)

    testthat::expect_equal(labeller(""), "")

    result <- labeller("___")
    testthat::expect_true(grepl("^\\s+$", result))

    result <- labeller("_start_end_")
    testthat::expect_true(grepl("start", result))
    testthat::expect_true(grepl("end", result))
    testthat::expect_false(grepl("_", result, fixed = TRUE))
})

testthat::test_that(".defaultLabeller with different n values", {
    labeller_5 <- .defaultLabeller(n = 5)
    labeller_20 <- .defaultLabeller(n = 20)

    str <- "test_string_here"
    result_5 <- labeller_5(str)
    result_20 <- labeller_20(str)

    testthat::expect_false(grepl("_", result_5, fixed = TRUE))
    testthat::expect_false(grepl("_", result_20, fixed = TRUE))
    testthat::expect_type(result_5, "character")
    testthat::expect_type(result_20, "character")
})
