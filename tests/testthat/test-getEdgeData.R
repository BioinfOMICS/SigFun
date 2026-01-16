library(testthat)
library(igraph)

test_that(".getEdgeData returns basic edge list", {
    g <- make_graph(~ A-B, B-C, C-A)

    result <- .getEdgeData(g)

    expect_s3_class(result, "data.frame")
    expect_equal(ncol(result), 2)
    expect_equal(nrow(result), 3)
})

test_that(".getEdgeData includes edge attributes", {
    g <- make_graph(~ A-B, B-C)
    E(g)$weight <- c(1.5, 2.5)
    E(g)$type <- c("strong", "weak")

    result <- .getEdgeData(g)

    expect_equal(ncol(result), 4)
    expect_true("weight" %in% names(result))
    expect_true("type" %in% names(result))
    expect_equal(result$weight, c(1.5, 2.5))
    expect_equal(result$type, c("strong", "weak"))
})

test_that(".getEdgeData handles graph with no edge attributes", {
    g <- make_graph(~ A-B, B-C, C-D)

    result <- .getEdgeData(g)

    expect_equal(ncol(result), 2)
    expect_equal(nrow(result), 3)
})

test_that(".getEdgeData handles empty graph", {
    g <- make_empty_graph(n = 3)

    result <- .getEdgeData(g)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 0)
    expect_equal(ncol(result), 2)
})

test_that(".getEdgeData handles single edge", {
    g <- make_graph(~ A-B)
    E(g)$label <- "connection"

    result <- .getEdgeData(g)

    expect_equal(nrow(result), 1)
    expect_equal(result$label, "connection")
})

test_that(".getEdgeData preserves multiple edge attributes", {
    g <- make_graph(~ A-B, B-C, C-D)
    E(g)$weight <- c(1, 2, 3)
    E(g)$color <- c("red", "blue", "green")
    E(g)$width <- c(0.5, 1.0, 1.5)

    result <- .getEdgeData(g)

    expect_equal(ncol(result), 5)
    expect_equal(result$weight, c(1, 2, 3))
    expect_equal(result$color, c("red", "blue", "green"))
    expect_equal(result$width, c(0.5, 1.0, 1.5))
})

test_that(".getEdgeData handles directed graph", {
    g <- make_graph(~ A-+B, B-+C, C-+A)
    E(g)$direction <- c("forward", "forward", "forward")

    result <- .getEdgeData(g)

    expect_equal(nrow(result), 3)
    expect_true("direction" %in% names(result))
})

test_that(".getEdgeData column names are correct", {
    g <- make_graph(~ X-Y, Y-Z)

    result <- .getEdgeData(g)

    expect_equal(names(result), c("V1", "V2"))
})

test_that(".getEdgeData handles NA values in edge attributes", {
    g <- make_graph(~ A-B, B-C)
    E(g)$weight <- c(1.0, NA)
    E(g)$label <- c(NA, "test")

    result <- .getEdgeData(g)

    expect_true(is.na(result$weight[2]))
    expect_true(is.na(result$label[1]))
    expect_equal(result$weight[1], 1.0)
    expect_equal(result$label[2], "test")
})

test_that(".getEdgeData handles complex graph structure", {
    g <- graph_from_data_frame(
        data.frame(from = c("A", "B", "C", "D"),
                   to = c("B", "C", "D", "A"))
    )
    E(g)$id <- 1:4

    result <- .getEdgeData(g)

    expect_equal(nrow(result), 4)
    expect_equal(result$id, 1:4)
})
