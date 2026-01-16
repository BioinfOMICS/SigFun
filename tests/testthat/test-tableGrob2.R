library(testthat)

test_that(".tableGrob2 returns gtable when p is NULL", {
    skip_if_not_installed("gridExtra")

    d <- data.frame(
        category = c("A", "B", "C"),
        value = c(10, 20, 30)
    )
    rownames(d) <- c("cat1", "cat2", "cat3")

    result <- .tableGrob2(d, p = NULL)

    expect_s3_class(result, "gtable")
})

test_that(".tableGrob2 sorts data by rownames", {
    skip_if_not_installed("gridExtra")

    d <- data.frame(value = c(10, 20, 30))
    rownames(d) <- c("C", "A", "B")

    result <- .tableGrob2(d, p = NULL)

    expect_s3_class(result, "gtable")
})

test_that(".tableGrob2 handles single row data frame", {
    skip_if_not_installed("gridExtra")

    d <- data.frame(value = 100)
    rownames(d) <- "single"

    result <- .tableGrob2(d, p = NULL)

    expect_s3_class(result, "gtable")
})

test_that(".tableGrob2 handles multiple columns", {
    skip_if_not_installed("gridExtra")

    d <- data.frame(
        col1 = c(1, 2, 3),
        col2 = c("A", "B", "C"),
        col3 = c(TRUE, FALSE, TRUE)
    )
    rownames(d) <- c("r1", "r2", "r3")

    result <- .tableGrob2(d, p = NULL)

    expect_s3_class(result, "gtable")
})

test_that(".tableGrob2 has expected grob structure", {
    skip_if_not_installed("gridExtra")

    d <- data.frame(x = 1:5, y = letters[1:5])
    rownames(d) <- paste0("row", 1:5)

    result <- .tableGrob2(d, p = NULL)

    expect_s3_class(result, "gtable")
    expect_true("layout" %in% names(result))
    expect_true("grobs" %in% names(result))
})

test_that(".tableGrob2 handles unsorted rownames", {
    skip_if_not_installed("gridExtra")

    d <- data.frame(value = c(30, 10, 20))
    rownames(d) <- c("Z", "A", "M")

    result <- .tableGrob2(d, p = NULL)

    expect_s3_class(result, "gtable")
})

test_that(".tableGrob2 works with numeric rownames", {
    skip_if_not_installed("gridExtra")

    d <- data.frame(data = c(100, 200, 300))
    rownames(d) <- c("3", "1", "2")

    result <- .tableGrob2(d, p = NULL)

    expect_s3_class(result, "gtable")
})

test_that(".tableGrob2 layout is a data frame", {
    skip_if_not_installed("gridExtra")

    d <- data.frame(a = 1:3, b = 4:6)
    rownames(d) <- c("x", "y", "z")

    result <- .tableGrob2(d, p = NULL)

    expect_true(is.data.frame(result$layout))
    expect_true("name" %in% names(result$layout))
})

test_that(".tableGrob2 handles empty data frame", {
    skip_if_not_installed("gridExtra")

    d <- data.frame()

    expect_error(
        .tableGrob2(d, p = NULL),
        "arguments imply differing number of rows"
    )
})
