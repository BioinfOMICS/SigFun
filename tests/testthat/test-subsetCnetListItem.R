library(testthat)

test_that(".subsetCnetListItem returns all when showItem is 'all'", {
    x <- list(
        gene1 = c("GO:0001", "GO:0002", "GO:0003"),
        gene2 = c("GO:0004", "GO:0005"),
        gene3 = c("GO:0001", "GO:0006")
    )

    result <- .subsetCnetListItem(x, "all")
    expect_equal(result, x)
    expect_equal(length(result), 3)
})

test_that(".subsetCnetListItem filters items correctly", {
    x <- list(
        gene1 = c("GO:0001", "GO:0002", "GO:0003"),
        gene2 = c("GO:0004", "GO:0005"),
        gene3 = c("GO:0001", "GO:0006")
    )

    result <- .subsetCnetListItem(x, "GO:0001")
    expect_equal(result$gene1, "GO:0001")
    expect_equal(result$gene2, character(0))
    expect_equal(result$gene3, "GO:0001")

    result <- .subsetCnetListItem(x, c("GO:0001", "GO:0005"))
    expect_equal(result$gene1, "GO:0001")
    expect_equal(result$gene2, "GO:0005")
    expect_equal(result$gene3, "GO:0001")
})

test_that(".subsetCnetListItem handles non-matching items", {
    x <- list(
        gene1 = c("GO:0001", "GO:0002"),
        gene2 = c("GO:0003", "GO:0004")
    )

    result <- .subsetCnetListItem(x, "GO:9999")
    expect_equal(result$gene1, character(0))
    expect_equal(result$gene2, character(0))

    result <- .subsetCnetListItem(x, c("GO:9999", "GO:8888"))
    expect_equal(result$gene1, character(0))
    expect_equal(result$gene2, character(0))
})

test_that(".subsetCnetListItem preserves list structure", {
    x <- list(
        gene1 = c("A", "B", "C"),
        gene2 = c("D", "E"),
        gene3 = c("F")
    )

    result <- .subsetCnetListItem(x, c("A", "D", "F"))
    expect_type(result, "list")
    expect_equal(length(result), 3)
    expect_equal(names(result), c("gene1", "gene2", "gene3"))
    expect_equal(result$gene1, "A")
    expect_equal(result$gene2, "D")
    expect_equal(result$gene3, "F")
})

test_that(".subsetCnetListItem handles edge cases", {
    x <- list()
    result <- .subsetCnetListItem(x, "all")
    expect_equal(result, list())

    x <- list(gene1 = character(0))
    result <- .subsetCnetListItem(x, "GO:0001")
    expect_equal(result$gene1, character(0))

    x <- list(gene1 = c("A", "B"))
    result <- .subsetCnetListItem(x, character(0))
    expect_equal(result$gene1, character(0))
})

test_that(".subsetCnetListItem handles partial matches", {
    x <- list(
        gene1 = c("GO:0001", "GO:0002", "GO:0003"),
        gene2 = c("GO:0002", "GO:0004"),
        gene3 = c("GO:0005")
    )

    result <- .subsetCnetListItem(x, c("GO:0002", "GO:0005"))
    expect_equal(result$gene1, "GO:0002")
    expect_equal(result$gene2, "GO:0002")
    expect_equal(result$gene3, "GO:0005")
})
