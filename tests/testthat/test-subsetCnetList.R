library(testthat)

test_that(".subsetCnetList filters by character vector", {
    x <- list(a = 1:3, b = 4:6, c = 7:9, d = 10:12)

    result <- .subsetCnetList(x, c("a", "c"))
    expect_equal(names(result), c("a", "c"))
    expect_equal(result$a, 1:3)
    expect_equal(result$c, 7:9)

    result <- .subsetCnetList(x, c("e", "f"))
    expect_equal(length(result), 0)

    result <- .subsetCnetList(x, c("a", "e", "d"))
    expect_equal(names(result), c("a", "d"))
})

test_that(".subsetCnetList filters by numeric index", {
    x <- list(a = 1:3, b = 4:6, c = 7:9, d = 10:12)

    result <- .subsetCnetList(x, 2)
    expect_equal(length(result), 2)
    expect_equal(names(result), c("a", "b"))

    expect_warning(
        result <- .subsetCnetList(x, c(1, 3, 4)),
        "coercing argument of type 'double' to logical"
    )
    expect_equal(length(result), 3)
    expect_equal(names(result), c("a", "c", "d"))

    result <- .subsetCnetList(x, integer(0))
    expect_equal(length(result), 0)
})

test_that(".subsetCnetList handles out of range indices", {
    x <- list(a = 1:3, b = 4:6, c = 7:9)

    expect_warning(
        result <- .subsetCnetList(x, c(1, 2, 5, 10)),
        "coercing argument of type 'double' to logical"
    )
    expect_equal(length(result), 2)
    expect_equal(names(result), c("a", "b"))

    expect_warning(
        result <- .subsetCnetList(x, c(5, 10)),
        "coercing argument of type 'double' to logical"
    )
    expect_equal(length(result), 0)
})

test_that(".subsetCnetList handles edge cases", {
    x <- list()
    result <- .subsetCnetList(x, character(0))
    expect_equal(length(result), 0)

    x <- list(a = 1:3)
    result <- .subsetCnetList(x, 1)
    expect_equal(length(result), 1)
    expect_equal(names(result), "a")
})

test_that(".subsetCnetList preserves data structure", {
    x <- list(
        gene1 = c("GO:0001", "GO:0002"),
        gene2 = c("GO:0003"),
        gene3 = c("GO:0001", "GO:0004", "GO:0005")
    )

    result <- .subsetCnetList(x, c("gene1", "gene3"))
    expect_type(result, "list")
    expect_equal(result$gene1, c("GO:0001", "GO:0002"))
    expect_equal(result$gene3, c("GO:0001", "GO:0004", "GO:0005"))
})
