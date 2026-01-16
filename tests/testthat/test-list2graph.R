library(testthat)
library(igraph)

test_that(".list2graph creates undirected graph by default", {
    inputList <- list(
        cat1 = c("gene1", "gene2"),
        cat2 = c("gene2", "gene3")
    )

    g <- .list2graph(inputList)

    expect_s3_class(g, "igraph")
    expect_false(is_directed(g))
    expect_equal(vcount(g), 5)
    expect_equal(ecount(g), 4)
})

test_that(".list2graph creates directed graph when specified", {
    inputList <- list(
        cat1 = c("gene1", "gene2"),
        cat2 = c("gene3")
    )

    g <- .list2graph(inputList, directed = TRUE)

    expect_s3_class(g, "igraph")
    expect_true(is_directed(g))
})

test_that(".list2graph sets .isCategory attribute correctly", {
    inputList <- list(
        cat1 = c("gene1", "gene2"),
        cat2 = c("gene3")
    )

    g <- .list2graph(inputList)

    expect_true(all(c(".isCategory") %in% vertex_attr_names(g)))

    category_vertices <- V(g)[V(g)$.isCategory]$name
    expect_true(all(c("cat1", "cat2") %in% category_vertices))

    item_vertices <- V(g)[!V(g)$.isCategory]$name
    expect_true(all(c("gene1", "gene2", "gene3") %in% item_vertices))
})

test_that(".list2graph sets size attribute based on degree", {
    inputList <- list(
        cat1 = c("gene1", "gene2"),
        cat2 = c("gene2")
    )

    g <- .list2graph(inputList)

    expect_true("size" %in% vertex_attr_names(g))

    gene2_idx <- which(V(g)$name == "gene2")
    expect_equal(as.numeric(V(g)$size[gene2_idx]), as.numeric(degree(g)[gene2_idx]))
})

test_that(".list2graph handles single category", {
    inputList <- list(
        cat1 = c("gene1", "gene2", "gene3")
    )

    g <- .list2graph(inputList)

    expect_equal(vcount(g), 4)
    expect_equal(ecount(g), 3)
    expect_equal(sum(V(g)$.isCategory), 1)
})

test_that(".list2graph handles empty list elements", {
    inputList <- list(
        cat1 = c("gene1"),
        cat2 = character(0),
        cat3 = c("gene2")
    )

    g <- .list2graph(inputList)

    expect_s3_class(g, "igraph")
    expect_true("cat1" %in% V(g)$name)
    expect_true("cat3" %in% V(g)$name)
})

test_that(".list2graph handles overlapping items", {
    inputList <- list(
        cat1 = c("gene1", "gene2"),
        cat2 = c("gene2", "gene3"),
        cat3 = c("gene1", "gene3")
    )

    g <- .list2graph(inputList)

    categories <- V(g)[V(g)$.isCategory]$name
    items <- V(g)[!V(g)$.isCategory]$name

    expect_equal(length(categories), 3)
    expect_equal(length(items), 3)
    expect_true(all(c("gene1", "gene2", "gene3") %in% items))
})

test_that(".list2graph vertex names are unique", {
    inputList <- list(
        cat1 = c("gene1", "gene2", "gene1"),
        cat2 = c("gene2", "gene3")
    )

    g <- .list2graph(inputList)

    vertex_names <- V(g)$name
    expect_equal(length(vertex_names), length(unique(vertex_names)))
})

test_that(".list2graph size attribute equals degree", {
    inputList <- list(
        cat1 = c("A", "B", "C"),
        cat2 = c("B", "D"),
        cat3 = c("C", "D", "E")
    )

    g <- .list2graph(inputList)

    expect_equal(as.numeric(V(g)$size), as.numeric(degree(g)))
})
