library(testthat)

test_that(".updateN returns character vector when showCategory is character", {
    x <- list(
        pathway1 = c("gene1", "gene2"),
        pathway2 = c("gene3", "gene4"),
        pathway3 = c("gene5", "gene6")
    )

    result <- .updateN(x, c("pathway1", "pathway3"))

    expect_type(result, "character")
    expect_equal(result, c("pathway1", "pathway3"))
})

test_that(".updateN filters non-existent names from character showCategory", {
    x <- list(
        setA = c("a", "b"),
        setB = c("c", "d")
    )

    result <- .updateN(x, c("setA", "setC", "setD"))

    expect_equal(result, "setA")
    expect_equal(length(result), 1)
})

test_that(".updateN works with data.frame and character showCategory", {
    x <- data.frame(
        ID = c("GO:001", "GO:002", "GO:003"),
        Description = c("Process 1", "Process 2", "Process 3")
    )

    result <- .updateN(x, c("Process 1", "Process 3"))

    expect_equal(result, c("Process 1", "Process 3"))
})

test_that(".updateN filters non-existent descriptions from data.frame", {
    x <- data.frame(
        ID = c("ID1", "ID2"),
        Description = c("Pathway A", "Pathway B")
    )

    result <- .updateN(x, c("Pathway A", "Pathway C", "Pathway D"))

    expect_equal(result, "Pathway A")
})

test_that(".updateN returns numeric when showCategory is numeric", {
    x <- list(a = 1, b = 2, c = 3, d = 4, e = 5)

    result <- .updateN(x, 3)

    expect_type(result, "double")
    expect_equal(result, 3)
})

test_that(".updateN limits n to list length when n is larger", {
    x <- list(item1 = "x", item2 = "y", item3 = "z")

    result <- .updateN(x, 10)

    expect_equal(result, 3)
})

test_that(".updateN limits n to data.frame rows when n is larger", {
    x <- data.frame(
        ID = c("A", "B"),
        Description = c("Desc A", "Desc B")
    )

    result <- .updateN(x, 5)

    expect_equal(result, 2)
})

test_that(".updateN returns n when n is smaller than list length", {
    x <- list(p1 = 1, p2 = 2, p3 = 3, p4 = 4)

    result <- .updateN(x, 2)

    expect_equal(result, 2)
})

test_that(".updateN returns n when n is smaller than data.frame rows", {
    x <- data.frame(
        ID = 1:10,
        Description = paste0("Item", 1:10)
    )

    result <- .updateN(x, 5)

    expect_equal(result, 5)
})

test_that(".updateN handles empty character vector", {
    x <- list(a = 1, b = 2)

    result <- .updateN(x, character(0))

    expect_equal(result, character(0))
    expect_equal(length(result), 0)
})

test_that(".updateN handles all non-matching character values", {
    x <- data.frame(Description = c("A", "B", "C"))

    result <- .updateN(x, c("X", "Y", "Z"))

    expect_equal(result, character(0))
})

test_that(".updateN handles n = 0 with list", {
    x <- list(a = 1, b = 2)

    result <- .updateN(x, 0)

    expect_equal(result, 0)
})

test_that(".updateN handles n = 1 with list", {
    x <- list(item1 = "x", item2 = "y")

    result <- .updateN(x, 1)

    expect_equal(result, 1)
})

test_that(".updateN handles single element list", {
    x <- list(only = "one")

    result <- .updateN(x, 5)

    expect_equal(result, 1)
})

test_that(".updateN handles single row data.frame", {
    x <- data.frame(ID = "A", Description = "Single")

    result <- .updateN(x, 10)

    expect_equal(result, 1)
})

test_that(".updateN preserves order of character showCategory", {
    x <- list(
        first = 1,
        second = 2,
        third = 3
    )

    result <- .updateN(x, c("third", "first"))

    expect_equal(result, c("third", "first"))
})
