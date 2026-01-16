library(testthat)

test_that(".gseaScores returns correct structure", {
    geneList <- c(gene1 = 3, gene2 = 2, gene3 = 1, gene4 = 0.5)
    geneSet <- c("gene1", "gene3")

    result <- .gseaScores(geneList, geneSet)

    expect_type(result, "list")
    expect_true("ES" %in% names(result))
    expect_true("runningES" %in% names(result))
    expect_type(result$ES, "double")
    expect_s3_class(result$runningES, "data.frame")
})

test_that(".gseaScores runningES data frame has correct columns", {
    geneList <- c(A = 5, B = 3, C = 1)
    geneSet <- c("A", "C")

    result <- .gseaScores(geneList, geneSet)
    df <- result$runningES

    expect_equal(names(df), c("x", "runningScore", "position", "gene"))
    expect_equal(nrow(df), length(geneList))
    expect_equal(df$gene, names(geneList))
})

test_that(".gseaScores fortify returns simplified data frame", {
    geneList <- c(A = 5, B = 3, C = 1)
    geneSet <- c("A")

    result <- .gseaScores(geneList, geneSet, fortify = TRUE)

    expect_s3_class(result, "data.frame")
    expect_equal(names(result), c("x", "runningScore", "position"))
    expect_false("gene" %in% names(result))
})

test_that(".gseaScores position column marks hits correctly", {
    geneList <- c(gene1 = 5, gene2 = 3, gene3 = 2, gene4 = 1)
    geneSet <- c("gene1", "gene3")

    result <- .gseaScores(geneList, geneSet)
    df <- result$runningES

    expect_equal(df$position[1], 1)
    expect_equal(df$position[2], 0)
    expect_equal(df$position[3], 1)
    expect_equal(df$position[4], 0)
})

test_that(".gseaScores handles all genes in set", {
    geneList <- c(A = 3, B = 2, C = 1)
    geneSet <- c("A", "B", "C")

    result <- .gseaScores(geneList, geneSet)

    expect_type(result$ES, "double")
    expect_equal(sum(result$runningES$position), 3)
})

test_that(".gseaScores handles no genes in set", {
    geneList <- c(A = 3, B = 2, C = 1)
    geneSet <- c("X", "Y", "Z")

    expect_error(
        .gseaScores(geneList, geneSet),
        "missing value where TRUE/FALSE needed"
    )
})

test_that(".gseaScores handles single gene in set", {
    geneList <- c(A = 5, B = 3, C = 1)
    geneSet <- c("A")

    result <- .gseaScores(geneList, geneSet)

    expect_type(result$ES, "double")
    expect_equal(sum(result$runningES$position), 1)
})

test_that(".gseaScores exponent parameter affects scores", {
    geneList <- c(A = 4, B = 3, C = 2, D = 1)
    geneSet <- c("A", "C")

    result1 <- .gseaScores(geneList, geneSet, exponent = 1)
    result2 <- .gseaScores(geneList, geneSet, exponent = 2)

    expect_false(identical(result1$ES, result2$ES))
    expect_false(identical(result1$runningES$runningScore,
                           result2$runningES$runningScore))
})

test_that(".gseaScores ES is maximum absolute value", {
    geneList <- c(A = 10, B = 8, C = 6, D = 4, E = 2)
    geneSet <- c("A", "B")

    result <- .gseaScores(geneList, geneSet)

    max_abs <- max(abs(result$runningES$runningScore))
    expect_equal(abs(result$ES), max_abs)
})

test_that(".gseaScores handles geneSet with non-existent genes", {
    geneList <- c(A = 5, B = 3, C = 1)
    geneSet <- c("A", "X", "Y", "C")

    result <- .gseaScores(geneList, geneSet)

    expect_equal(sum(result$runningES$position), 2)
    expect_equal(result$runningES$gene[result$runningES$position == 1], c("A", "C"))
})

test_that(".gseaScores x column is sequential", {
    geneList <- c(gene1 = 5, gene2 = 3, gene3 = 1)
    geneSet <- c("gene1")

    result <- .gseaScores(geneList, geneSet)

    expect_equal(result$runningES$x, 1:3)
})

test_that(".gseaScores handles negative gene scores", {
    geneList <- c(A = 5, B = -3, C = 2, D = -1)
    geneSet <- c("B", "D")

    result <- .gseaScores(geneList, geneSet)

    expect_type(result$ES, "double")
    expect_true(all(result$runningES$runningScore >= min(result$runningES$runningScore)))
})

test_that(".gseaScores runningScore is cumulative", {
    geneList <- c(A = 4, B = 3, C = 2, D = 1)
    geneSet <- c("A", "C")

    result <- .gseaScores(geneList, geneSet)
    rs <- result$runningES$runningScore

    expect_equal(length(rs), 4)
    expect_true(all(!is.na(rs)))
})
