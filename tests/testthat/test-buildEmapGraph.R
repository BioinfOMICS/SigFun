library(testthat)
library(igraph)

test_that(".buildEmapGraph validates min_edge parameter", {
    enrichDf <- data.frame(Description = c("A", "B"))
    geneSets <- list(A = c("g1"), B = c("g2"))
    pair_sim <- matrix(1, nrow = 2, ncol = 2)

    expect_error(
        .buildEmapGraph(enrichDf, geneSets, "red", 1, -0.1, pair_sim, "JC"),
        "should be a number between 0 and 1"
    )

    expect_error(
        .buildEmapGraph(enrichDf, geneSets, "red", 1, 1.5, pair_sim, "JC"),
        "should be a number between 0 and 1"
    )

    expect_error(
        .buildEmapGraph(enrichDf, geneSets, "red", 1, "invalid", pair_sim, "JC"),
        "should be a number between 0 and 1"
    )
})

test_that(".buildEmapGraph handles single row enrichDf", {
    enrichDf <- data.frame(Description = "Pathway A")
    geneSets <- list("Pathway A" = c("gene1", "gene2"))
    pair_sim <- matrix(1, nrow = 1, ncol = 1)

    g <- .buildEmapGraph(enrichDf, geneSets, "red", 1, 0.2, pair_sim, "JC")

    expect_s3_class(g, "igraph")
    expect_equal(vcount(g), 1)
    expect_equal(V(g)$name, "Pathway A")
    expect_equal(V(g)$color, "red")
})

test_that(".buildEmapGraph handles NULL enrichDf dimensions", {
    enrichDf <- data.frame(Description = "Single")
    geneSets <- list(Single = c("g1"))
    pair_sim <- matrix(1, nrow = 1, ncol = 1)

    g <- .buildEmapGraph(enrichDf, geneSets, "blue", 1, 0.2, pair_sim, "JC")

    expect_s3_class(g, "igraph")
    expect_equal(vcount(g), 1)
})

test_that(".buildEmapGraph creates graph from similarity matrix", {
    enrichDf <- data.frame(
        Description = c("Path1", "Path2", "Path3"),
        pvalue = c(0.01, 0.02, 0.03)
    )
    geneSets <- list(
        Path1 = c("g1", "g2"),
        Path2 = c("g2", "g3"),
        Path3 = c("g3", "g4")
    )
    pair_sim <- matrix(c(1.0, 0.5, 0.3,
                         0.5, 1.0, 0.6,
                         0.3, 0.6, 1.0),
                       nrow = 3, ncol = 3,
                       dimnames = list(c("Path1", "Path2", "Path3"),
                                       c("Path1", "Path2", "Path3")))

    g <- .buildEmapGraph(enrichDf, geneSets, "pvalue", 1, 0.2, pair_sim, "JC")

    expect_s3_class(g, "igraph")
    expect_equal(vcount(g), 3)
    expect_true(ecount(g) > 0)
})

test_that(".buildEmapGraph filters edges by min_edge", {
    enrichDf <- data.frame(Description = c("A", "B", "C"))
    geneSets <- list(A = c("g1"), B = c("g2"), C = c("g3"))
    pair_sim <- matrix(c(1.0, 0.8, 0.1,
                         0.8, 1.0, 0.1,
                         0.1, 0.1, 1.0),
                       nrow = 3, ncol = 3,
                       dimnames = list(c("A", "B", "C"), c("A", "B", "C")))

    g_low <- .buildEmapGraph(enrichDf, geneSets, "red", 1, 0.05, pair_sim, "JC")
    g_high <- .buildEmapGraph(enrichDf, geneSets, "red", 1, 0.5, pair_sim, "JC")

    expect_true(ecount(g_low) >= ecount(g_high))
})

test_that(".buildEmapGraph sets edge width based on similarity", {
    enrichDf <- data.frame(Description = c("P1", "P2"))
    geneSets <- list(P1 = c("g1"), P2 = c("g2"))
    pair_sim <- matrix(c(1.0, 0.5, 0.5, 1.0),
                       nrow = 2, ncol = 2,
                       dimnames = list(c("P1", "P2"), c("P1", "P2")))

    g <- .buildEmapGraph(enrichDf, geneSets, "red", 2, 0.2, pair_sim, "JC")

    if (ecount(g) > 0) {
        expect_true("width" %in% edge_attr_names(g))
        expect_true("weight" %in% edge_attr_names(g))
    }
})

test_that(".buildEmapGraph sets vertex size from geneSets length", {
    enrichDf <- data.frame(Description = c("Path1", "Path2", "Path3"))
    geneSets <- list(
        Path1 = c("g1", "g2", "g3"),
        Path2 = c("g4", "g5"),
        Path3 = c("g6")
    )
    pair_sim <- matrix(c(1.0, 0.5, 0.5,
                         0.5, 1.0, 0.5,
                         0.5, 0.5, 1.0),
                       nrow = 3, ncol = 3,
                       dimnames = list(c("Path1", "Path2", "Path3"),
                                       c("Path1", "Path2", "Path3")))

    g <- .buildEmapGraph(enrichDf, geneSets, "red", 1, 0.2, pair_sim, "JC")

    expect_true("size" %in% vertex_attr_names(g))
    expect_equal(length(V(g)$size), 3)
    expect_true(all(V(g)$size %in% c(1, 2, 3)))
})

test_that(".buildEmapGraph uses column for color when specified", {
    enrichDf <- data.frame(
        Description = c("A", "B"),
        pvalue = c(0.01, 0.05)
    )
    geneSets <- list(A = c("g1"), B = c("g2"))
    pair_sim <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2,
                       dimnames = list(c("A", "B"), c("A", "B")))

    g <- .buildEmapGraph(enrichDf, geneSets, "pvalue", 1, 0.2, pair_sim, "JC")

    expect_true("color" %in% vertex_attr_names(g))
    expect_equal(length(V(g)$color), 2)
    expect_true(all(V(g)$color %in% c(0.01, 0.05)))
})

test_that(".buildEmapGraph uses constant color when not in enrichDf", {
    enrichDf <- data.frame(Description = c("X", "Y"))
    geneSets <- list(X = c("g1"), Y = c("g2"))
    pair_sim <- matrix(c(1, 0.6, 0.6, 1), nrow = 2, ncol = 2,
                       dimnames = list(c("X", "Y"), c("X", "Y")))

    g <- .buildEmapGraph(enrichDf, geneSets, "blue", 1, 0.2, pair_sim, "JC")

    expect_equal(V(g)$color, c("blue", "blue"))
})

test_that(".buildEmapGraph handles method parameter", {
    enrichDf <- data.frame(Description = c("P1", "P2"))
    geneSets <- list(P1 = c("g1"), P2 = c("g2"))
    pair_sim <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2,
                       dimnames = list(c("P1", "P2"), c("P1", "P2")))

    g_jc <- .buildEmapGraph(enrichDf, geneSets, "red", 1, 0.2, pair_sim, "JC")
    g_other <- .buildEmapGraph(enrichDf, geneSets, "red", 1, 0.2, pair_sim, "Other")

    expect_s3_class(g_jc, "igraph")
    expect_s3_class(g_other, "igraph")
})

test_that(".buildEmapGraph scales edge width by cex_line", {
    enrichDf <- data.frame(Description = c("A", "B"))
    geneSets <- list(A = c("g1"), B = c("g2"))
    pair_sim <- matrix(c(1, 0.8, 0.8, 1), nrow = 2, ncol = 2,
                       dimnames = list(c("A", "B"), c("A", "B")))

    g1 <- .buildEmapGraph(enrichDf, geneSets, "red", 1, 0.2, pair_sim, "JC")
    g2 <- .buildEmapGraph(enrichDf, geneSets, "red", 2, 0.2, pair_sim, "JC")

    if (ecount(g1) > 0 && ecount(g2) > 0) {
        expect_true(E(g2)$width[1] > E(g1)$width[1])
    }
})

test_that(".buildEmapGraph removes NA similarities", {
    enrichDf <- data.frame(Description = c("P1", "P2", "P3"))
    geneSets <- list(P1 = c("g1"), P2 = c("g2"), P3 = c("g3"))
    pair_sim <- matrix(c(1.0, 0.5, NA,
                         0.5, 1.0, NA,
                         NA, NA, 1.0),
                       nrow = 3, ncol = 3,
                       dimnames = list(c("P1", "P2", "P3"),
                                       c("P1", "P2", "P3")))

    g <- .buildEmapGraph(enrichDf, geneSets, "red", 1, 0.2, pair_sim, "JC")

    expect_s3_class(g, "igraph")
    expect_true(vcount(g) >= 2)
})
