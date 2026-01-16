library(testthat)

test_that(".getGsdata returns single pathway data", {
    geneList <- c(gene1 = 5, gene2 = 3, gene3 = 2, gene4 = 1)

    result_df <- data.frame(
        ID = "pathway1",
        Description = "Test Pathway",
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- "pathway1"

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(pathway1 = c("gene1", "gene3")),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = result_df)

    result <- .getGsdata(object, "pathway1")

    expect_s3_class(result, "data.frame")
    expect_true(all(result$Description == "Test Pathway"))
    expect_equal(nrow(result), length(geneList))
})

test_that(".getGsdata combines multiple pathway data", {
    geneList <- c(A = 5, B = 4, C = 3, D = 2)

    result_df <- data.frame(
        ID = c("pathway1", "pathway2"),
        Description = c("Pathway 1", "Pathway 2"),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("pathway1", "pathway2")

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(
                      pathway1 = c("A", "C"),
                      pathway2 = c("B", "D")
                  ),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = result_df)

    result <- .getGsdata(object, c("pathway1", "pathway2"))

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), length(geneList) * 2)
    expect_true("Pathway 1" %in% result$Description)
    expect_true("Pathway 2" %in% result$Description)
})

test_that(".getGsdata preserves all columns", {
    geneList <- c(gene1 = 3, gene2 = 2, gene3 = 1)

    result_df <- data.frame(
        ID = c("pathway1", "pathway2"),
        Description = c("Path 1", "Path 2"),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("pathway1", "pathway2")

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(
                      pathway1 = c("gene1"),
                      pathway2 = c("gene2")
                  ),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = result_df)

    result <- .getGsdata(object, c("pathway1", "pathway2"))

    expected_cols <- c("x", "runningScore", "position", "ymin", "ymax",
                       "geneList", "gene", "Description")
    expect_true(all(expected_cols %in% names(result)))
})

test_that(".getGsdata handles single pathway with length 1 vector", {
    geneList <- c(A = 5, B = 3)

    result_df <- data.frame(
        ID = "pathway1",
        Description = "Single",
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- "pathway1"

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(pathway1 = c("A")),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = result_df)

    result <- .getGsdata(object, "pathway1")

    expect_equal(nrow(result), 2)
    expect_equal(unique(result$Description), "Single")
})

test_that(".getGsdata rows equal geneList length times pathway count", {
    geneList <- c(A = 4, B = 3, C = 2)

    result_df <- data.frame(
        ID = c("p1", "p2", "p3"),
        Description = c("Path 1", "Path 2", "Path 3"),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("p1", "p2", "p3")

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(
                      p1 = c("A"),
                      p2 = c("B"),
                      p3 = c("C")
                  ),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = result_df)

    result <- .getGsdata(object, c("p1", "p2", "p3"))

    expect_equal(nrow(result), length(geneList) * 3)
})

test_that(".getGsdata handles numeric geneSetID", {
    geneList <- c(A = 3, B = 2, C = 1)

    result_df <- data.frame(
        ID = c("pathway1", "pathway2"),
        Description = c("First", "Second"),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("pathway1", "pathway2")

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(
                      pathway1 = c("A"),
                      pathway2 = c("B")
                  ),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = result_df)

    result_single <- .getGsdata(object, 1)
    result_multi <- .getGsdata(object, c(1, 2))

    expect_equal(unique(result_single$Description), "First")
    expect_equal(nrow(result_multi), length(geneList) * 2)
})

test_that(".getGsdata maintains data consistency across pathways", {
    geneList <- c(gene1 = 5, gene2 = 3)

    result_df <- data.frame(
        ID = c("p1", "p2"),
        Description = c("Desc1", "Desc2"),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("p1", "p2")

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(
                      p1 = c("gene1"),
                      p2 = c("gene2")
                  ),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = result_df)

    result <- .getGsdata(object, c("p1", "p2"))

    p1_rows <- result[result$Description == "Desc1", ]
    p2_rows <- result[result$Description == "Desc2", ]

    expect_equal(nrow(p1_rows), 2)
    expect_equal(nrow(p2_rows), 2)
    expect_equal(as.numeric(p1_rows$geneList), as.numeric(geneList))
    expect_equal(as.numeric(p2_rows$geneList), as.numeric(geneList))
})

test_that(".getGsdata works with gene2Symbol mapping", {
    geneList <- c(ENSG001 = 5, ENSG002 = 3, ENSG003 = 1)
    gene2Symbol <- c(ENSG001 = "GeneA", ENSG002 = "GeneB", ENSG003 = "GeneC")

    result_df <- data.frame(
        ID = c("pathway1", "pathway2"),
        Description = c("P1", "P2"),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("pathway1", "pathway2")

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(
                      pathway1 = c("ENSG001"),
                      pathway2 = c("ENSG002")
                  ),
                  params = list(exponent = 1),
                  gene2Symbol = gene2Symbol,
                  result = result_df)

    result <- .getGsdata(object, c("pathway1", "pathway2"))

    expect_true(all(c("GeneA", "GeneB", "GeneC") %in% result$gene))
})
