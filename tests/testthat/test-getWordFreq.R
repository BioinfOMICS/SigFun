library(testthat)

testthat::test_that(".getWordFreq counts word frequencies correctly", {
    wordd <- c("apple banana", "banana cherry", "apple banana cherry")
    result <- .getWordFreq(wordd)

    testthat::expect_type(result, "double")
    testthat::expect_true(all(result > 0))
    testthat::expect_equal(length(result), 3)
    testthat::expect_true(all(result[-1] <= result[-length(result)]))
})

testthat::test_that(".getWordFreq handles single word", {
    wordd <- c("apple", "apple", "apple")
    result <- .getWordFreq(wordd)

    testthat::expect_equal(length(result), 1)
    testthat::expect_equal(names(result), "apple")
    testthat::expect_equal(as.vector(result), 3)
})

testthat::test_that(".getWordFreq handles multiple occurrences", {
    wordd <- c("a b", "b c", "a b c")
    result <- .getWordFreq(wordd)

    testthat::expect_equal(result["b"], c(b = 3))
    testthat::expect_equal(result["a"], c(a = 2))
    testthat::expect_equal(result["c"], c(c = 2))
})

testthat::test_that(".getWordFreq handles empty input", {
    wordd <- character(0)
    result <- .getWordFreq(wordd)

    testthat::expect_equal(length(result), 0)
})

testthat::test_that(".getWordFreq handles single element", {
    wordd <- "word1 word2 word1"
    result <- .getWordFreq(wordd)

    testthat::expect_equal(result["word1"], c(word1 = 1))
    testthat::expect_equal(result["word2"], c(word2 = 1))
})

testthat::test_that(".getWordFreq orders by document frequency", {
    wordd <- c("a b c", "b c d", "c d e")
    result <- .getWordFreq(wordd)

    testthat::expect_equal(names(result)[1], "c")
    testthat::expect_equal(result[1], c(c = 3))
})
