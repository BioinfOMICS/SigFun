library(testthat)

test_that(".gsInfo returns correct data frame structure", {
    geneList <- c(gene1 = 5, gene2 = 3, gene3 = 2, gene4 = 1)

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(pathway1 = c("gene1", "gene3")),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = data.frame(
                      ID = "pathway1",
                      Description = "Test Pathway",
                      stringsAsFactors = FALSE
                  ))

    result <- .gsInfo(object, "pathway1")

    expect_s3_class(result, "data.frame")
    expect_true(all(c("x", "runningScore", "position", "ymin", "ymax",
                      "geneList", "gene", "Description") %in% names(result)))
})

test_that(".gsInfo accepts numeric geneSetID", {
    geneList <- c(A = 4, B = 3, C = 2)

    result_df <- data.frame(
        ID = c("pathway1", "pathway2"),
        Description = c("Path 1", "Path 2"),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("pathway1", "pathway2")

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(pathway1 = c("A"), pathway2 = c("B")),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = result_df)

    result <- .gsInfo(object, 1)

    expect_equal(unique(result$Description), "Path 1")
})

test_that(".gsInfo accepts character geneSetID", {
    geneList <- c(A = 4, B = 3, C = 2)

    result_df <- data.frame(
        ID = c("pathway1", "pathway2"),
        Description = c("Path 1", "Path 2"),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("pathway1", "pathway2")

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(pathway1 = c("A"), pathway2 = c("B")),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = result_df)

    result <- .gsInfo(object, "pathway2")

    expect_equal(unique(result$Description), "Path 2")
})

test_that(".gsInfo sets ymin and ymax for hits", {
    geneList <- c(gene1 = 5, gene2 = 3, gene3 = 1)

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(pathway1 = c("gene1", "gene3")),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = data.frame(
                      ID = "pathway1",
                      Description = "Test",
                      stringsAsFactors = FALSE
                  ))

    result <- .gsInfo(object, "pathway1")

    hits <- which(result$position == 1)
    non_hits <- which(result$position == 0)

    expect_true(all(result$ymin[hits] < 0))
    expect_true(all(result$ymax[hits] > 0))
    expect_true(all(result$ymin[non_hits] == 0))
    expect_true(all(result$ymax[non_hits] == 0))
})

test_that(".gsInfo calculates ymin/ymax based on runningScore range", {
    geneList <- c(A = 10, B = 8, C = 6, D = 4)

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(pathway1 = c("A")),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = data.frame(
                      ID = "pathway1",
                      Description = "Test",
                      stringsAsFactors = FALSE
                  ))

    result <- .gsInfo(object, "pathway1")

    h <- diff(range(result$runningScore)) / 20
    hits <- which(result$position == 1)

    expect_equal(result$ymin[hits], rep(-h, length(hits)))
    expect_equal(result$ymax[hits], rep(h, length(hits)))
})

test_that(".gsInfo includes geneList column", {
    geneList <- c(gene1 = 5, gene2 = 3, gene3 = 1)

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(pathway1 = c("gene1")),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = data.frame(
                      ID = "pathway1",
                      Description = "Test",
                      stringsAsFactors = FALSE
                  ))

    result <- .gsInfo(object, "pathway1")

    expect_equal(as.numeric(result$geneList), as.numeric(geneList))
})

test_that(".gsInfo uses gene names when gene2Symbol is empty", {
    geneList <- c(GENE1 = 5, GENE2 = 3, GENE3 = 1)

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(pathway1 = c("GENE1")),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = data.frame(
                      ID = "pathway1",
                      Description = "Test",
                      stringsAsFactors = FALSE
                  ))

    result <- .gsInfo(object, "pathway1")

    expect_equal(result$gene, names(geneList))
})

test_that(".gsInfo uses gene2Symbol mapping when available", {
    geneList <- c(ENSG001 = 5, ENSG002 = 3, ENSG003 = 1)
    gene2Symbol <- c(ENSG001 = "GeneA", ENSG002 = "GeneB", ENSG003 = "GeneC")

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(pathway1 = c("ENSG001")),
                  params = list(exponent = 1),
                  gene2Symbol = gene2Symbol,
                  result = data.frame(
                      ID = "pathway1",
                      Description = "Test",
                      stringsAsFactors = FALSE
                  ))

    result <- .gsInfo(object, "pathway1")

    expect_equal(result$gene, c("GeneA", "GeneB", "GeneC"))
})

test_that(".gsInfo adds Description column", {
    geneList <- c(A = 5, B = 3)

    result_df <- data.frame(
        ID = "pathway1",
        Description = "Important Pathway",
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- "pathway1"

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(pathway1 = c("A")),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = result_df)

    result <- .gsInfo(object, "pathway1")

    expect_true(all(result$Description == "Important Pathway"))
})

test_that(".gsInfo respects exponent parameter", {
    geneList <- c(A = 4, B = 3, C = 2)

    object1 <- new("gseaResult",
                   geneList = geneList,
                   geneSets = list(pathway1 = c("A", "C")),
                   params = list(exponent = 1),
                   gene2Symbol = character(0),
                   result = data.frame(
                       ID = "pathway1",
                       Description = "Test",
                       stringsAsFactors = FALSE
                   ))

    object2 <- new("gseaResult",
                   geneList = geneList,
                   geneSets = list(pathway1 = c("A", "C")),
                   params = list(exponent = 2),
                   gene2Symbol = character(0),
                   result = data.frame(
                       ID = "pathway1",
                       Description = "Test",
                       stringsAsFactors = FALSE
                   ))

    result1 <- .gsInfo(object1, "pathway1")
    result2 <- .gsInfo(object2, "pathway1")

    expect_false(identical(result1$runningScore, result2$runningScore))
})

test_that(".gsInfo handles single gene in pathway", {
    geneList <- c(A = 5, B = 3, C = 1)

    object <- new("gseaResult",
                  geneList = geneList,
                  geneSets = list(pathway1 = c("A")),
                  params = list(exponent = 1),
                  gene2Symbol = character(0),
                  result = data.frame(
                      ID = "pathway1",
                      Description = "Single Gene Pathway",
                      stringsAsFactors = FALSE
                  ))

    result <- .gsInfo(object, "pathway1")

    expect_equal(sum(result$position), 1)
    expect_equal(nrow(result), 3)
})
