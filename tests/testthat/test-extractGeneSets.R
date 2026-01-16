library(testthat)

test_that(".extractGeneSets handles list input", {
    geneSets <- list(
        pathway1 = c("gene1", "gene2", "gene3"),
        pathway2 = c("gene4", "gene5"),
        pathway3 = c("gene6", "gene7", "gene8")
    )

    result <- .extractGeneSets(geneSets, 2)

    expect_type(result, "list")
    expect_equal(length(result), 2)
    expect_equal(names(result), c("pathway1", "pathway2"))
})

test_that(".extractGeneSets handles numeric n with list", {
    geneSets <- list(
        set1 = c("A", "B"),
        set2 = c("C", "D"),
        set3 = c("E", "F"),
        set4 = c("G", "H")
    )

    result <- .extractGeneSets(geneSets, 3)

    expect_equal(length(result), 3)
    expect_equal(names(result), c("set1", "set2", "set3"))
})

test_that(".extractGeneSets handles character n with list", {
    geneSets <- list(
        pathway1 = c("gene1", "gene2"),
        pathway2 = c("gene3", "gene4"),
        pathway3 = c("gene5", "gene6")
    )

    result <- .extractGeneSets(geneSets, c("pathway1", "pathway3"))

    expect_equal(length(result), 2)
    expect_equal(names(result), c("pathway1", "pathway3"))
})

test_that(".extractGeneSets handles single character n", {
    geneSets <- list(
        setA = c("x", "y"),
        setB = c("z", "w")
    )

    result <- .extractGeneSets(geneSets, "setB")

    expect_equal(length(result), 1)
    expect_equal(names(result), "setB")
    expect_equal(result$setB, c("z", "w"))
})

test_that(".extractGeneSets handles n = 1", {
    geneSets <- list(
        first = c("a", "b", "c"),
        second = c("d", "e")
    )

    result <- .extractGeneSets(geneSets, 1)

    expect_equal(length(result), 1)
    expect_equal(names(result), "first")
})

test_that(".extractGeneSets handles all elements", {
    geneSets <- list(
        p1 = c("g1"),
        p2 = c("g2"),
        p3 = c("g3")
    )

    result <- .extractGeneSets(geneSets, 3)

    expect_equal(length(result), 3)
    expect_equal(result, geneSets)
})

test_that(".extractGeneSets preserves list structure", {
    geneSets <- list(
        pathway1 = c("gene1", "gene2", "gene3"),
        pathway2 = c("gene4"),
        pathway3 = c("gene5", "gene6")
    )

    result <- .extractGeneSets(geneSets, 2)

    expect_type(result, "list")
    expect_equal(result$pathway1, c("gene1", "gene2", "gene3"))
    expect_equal(result$pathway2, c("gene4"))
})

test_that(".extractGeneSets handles empty selection", {
    geneSets <- list(
        set1 = c("a", "b"),
        set2 = c("c", "d")
    )

    result <- .extractGeneSets(geneSets, character(0))

    expect_equal(length(result), 0)
})

test_that(".extractGeneSets handles non-existent names", {
    geneSets <- list(
        pathway1 = c("gene1"),
        pathway2 = c("gene2")
    )

    result <- .extractGeneSets(geneSets, c("pathway1", "nonexistent"))

    expect_equal(length(result), 1)
    expect_equal(names(result), "pathway1")
})

test_that(".extractGeneSets with enrichResult object extracts correctly", {
    skip_if_not_installed("DOSE")

    result_df <- data.frame(
        ID = c("GO:0001", "GO:0002", "GO:0003"),
        Description = c("Process 1", "Process 2", "Process 3"),
        GeneRatio = c("2/100", "2/100", "2/100"),
        BgRatio = c("10/1000", "10/1000", "10/1000"),
        pvalue = c(0.01, 0.02, 0.03),
        p.adjust = c(0.01, 0.02, 0.03),
        qvalue = c(0.01, 0.02, 0.03),
        geneID = c("gene1/gene2", "gene3/gene4", "gene5/gene6"),
        Count = c(2, 2, 2),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("GO:0001", "GO:0002", "GO:0003")

    mock_result <- new("enrichResult",
                       result = result_df,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,
                       organism = "test",
                       ontology = "test",
                       gene = c("gene1", "gene2", "gene3", "gene4", "gene5", "gene6"),
                       keytype = "SYMBOL",
                       universe = character(0),
                       geneSets = list(
                           "GO:0001" = c("gene1", "gene2"),
                           "GO:0002" = c("gene3", "gene4"),
                           "GO:0003" = c("gene5", "gene6")
                       ),
                       readable = FALSE)

    result <- .extractGeneSets(mock_result, 2)

    expect_type(result, "list")
    expect_equal(length(result), 2)
    expect_equal(names(result), c("Process 1", "Process 2"))
})

test_that(".extractGeneSets with enrichResult uses Description as names", {
    skip_if_not_installed("DOSE")

    result_df <- data.frame(
        ID = c("ID001", "ID002"),
        Description = c("Pathway A", "Pathway B"),
        GeneRatio = c("2/100", "2/100"),
        BgRatio = c("10/1000", "10/1000"),
        pvalue = c(0.01, 0.02),
        p.adjust = c(0.01, 0.02),
        qvalue = c(0.01, 0.02),
        geneID = c("g1/g2", "g3/g4"),
        Count = c(2, 2),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("ID001", "ID002")

    mock_result <- new("enrichResult",
                       result = result_df,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,
                       organism = "test",
                       ontology = "test",
                       gene = c("g1", "g2", "g3", "g4"),
                       keytype = "SYMBOL",
                       universe = character(0),
                       geneSets = list(
                           ID001 = c("g1", "g2"),
                           ID002 = c("g3", "g4")
                       ),
                       readable = FALSE)

    result <- .extractGeneSets(mock_result, 1)

    expect_equal(names(result), "Pathway A")
    expect_equal(result[[1]], c("g1", "g2"))
})
