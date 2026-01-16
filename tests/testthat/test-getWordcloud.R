library(testthat)

testthat::test_that(".getWordcloud extracts top words from cluster", {
    ggData <- data.frame(
        name = c("gene A in pathway 1", "gene B in pathway 2",
                 "gene C in pathway 1", "gene D in pathway 2"),
        color2 = c("cluster1", "cluster1", "cluster2", "cluster2"),
        stringsAsFactors = FALSE
    )

    result <- .getWordcloud(cluster = "cluster1", ggData = ggData, nWords = 2)

    testthat::expect_type(result, "character")
    testthat::expect_true(nchar(result) > 0)
    testthat::expect_true(grepl("gene|pathway", result))
})

testthat::test_that(".getWordcloud respects nWords limit", {
    ggData <- data.frame(
        name = c("word1 word2 word3", "word1 word2", "word4 word5"),
        color2 = c("A", "A", "B"),
        stringsAsFactors = FALSE
    )

    result <- .getWordcloud(cluster = "A", ggData = ggData, nWords = 1)
    words <- strsplit(result, " ")[[1]]

    testthat::expect_lte(length(words), 1)
})

testthat::test_that(".getWordcloud handles single entry cluster", {
    ggData <- data.frame(
        name = c("test pathway", "other pathway"),
        color2 = c("X", "Y"),
        stringsAsFactors = FALSE
    )

    result <- .getWordcloud(cluster = "X", ggData = ggData, nWords = 5)

    testthat::expect_type(result, "character")
})

testthat::test_that(".getWordcloud removes common words", {
    ggData <- data.frame(
        name = c("gene 1 in pathway and test", "gene 2 of system"),
        color2 = c("A", "A"),
        stringsAsFactors = FALSE
    )

    result <- .getWordcloud(cluster = "A", ggData = ggData, nWords = 5)

    testthat::expect_false(grepl(" in ", result, fixed = TRUE))
    testthat::expect_false(grepl(" and ", result, fixed = TRUE))
    testthat::expect_false(grepl(" of ", result, fixed = TRUE))
})

testthat::test_that(".getWordcloud handles multiple clusters", {
    ggData <- data.frame(
        name = rep(c("alpha beta", "gamma delta", "epsilon zeta"), 2),
        color2 = rep(c("C1", "C2", "C3"), 2),
        stringsAsFactors = FALSE
    )

    result1 <- .getWordcloud(cluster = "C1", ggData = ggData, nWords = 2)
    result2 <- .getWordcloud(cluster = "C2", ggData = ggData, nWords = 2)

    testthat::expect_type(result1, "character")
    testthat::expect_type(result2, "character")
})

testthat::test_that(".getWordcloud orders words by rank", {
    ggData <- data.frame(
        name = c("word1 word2 word3", "word1 word2", "word1"),
        color2 = c("G", "G", "G"),
        stringsAsFactors = FALSE
    )

    result <- .getWordcloud(cluster = "G", ggData = ggData, nWords = 3)

    testthat::expect_type(result, "character")
    testthat::expect_true(grepl("word1", result))
})
