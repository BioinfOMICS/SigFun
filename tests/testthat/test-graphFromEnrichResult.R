library(testthat)
library(igraph)

test_that(".graphFromEnrichResult returns list with graph and geneSet", {
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

    termsim <- matrix(c(1.0, 0.5, 0.3,
                        0.5, 1.0, 0.4,
                        0.3, 0.4, 1.0),
                      nrow = 3, ncol = 3,
                      dimnames = list(c("Process 1", "Process 2", "Process 3"),
                                      c("Process 1", "Process 2", "Process 3")))

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
                       readable = FALSE,
                       termsim = termsim,
                       method = "JC")

    result <- .graphFromEnrichResult(mock_result, 2)

    expect_type(result, "list")
    expect_true("graph" %in% names(result))
    expect_true("geneSet" %in% names(result))
    expect_s3_class(result$graph, "igraph")
    expect_type(result$geneSet, "list")
})

test_that(".graphFromEnrichResult handles numeric showCategory", {
    skip_if_not_installed("DOSE")

    result_df <- data.frame(
        ID = c("ID1", "ID2", "ID3", "ID4"),
        Description = c("Path1", "Path2", "Path3", "Path4"),
        GeneRatio = c("2/100", "2/100", "2/100", "2/100"),
        BgRatio = c("10/1000", "10/1000", "10/1000", "10/1000"),
        pvalue = c(0.01, 0.02, 0.03, 0.04),
        p.adjust = c(0.01, 0.02, 0.03, 0.04),
        qvalue = c(0.01, 0.02, 0.03, 0.04),
        geneID = c("g1/g2", "g3/g4", "g5/g6", "g7/g8"),
        Count = c(2, 2, 2, 2),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("ID1", "ID2", "ID3", "ID4")

    termsim <- matrix(1, nrow = 4, ncol = 4,
                      dimnames = list(c("Path1", "Path2", "Path3", "Path4"),
                                      c("Path1", "Path2", "Path3", "Path4")))
    diag(termsim) <- 1

    mock_result <- new("enrichResult",
                       result = result_df,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,
                       organism = "test",
                       ontology = "test",
                       gene = paste0("g", 1:8),
                       keytype = "SYMBOL",
                       universe = character(0),
                       geneSets = list(
                           ID1 = c("g1", "g2"),
                           ID2 = c("g3", "g4"),
                           ID3 = c("g5", "g6"),
                           ID4 = c("g7", "g8")
                       ),
                       readable = FALSE,
                       termsim = termsim,
                       method = "JC")

    result <- .graphFromEnrichResult(mock_result, 2)

    expect_equal(length(result$geneSet), 2)
})

test_that(".graphFromEnrichResult handles character showCategory", {
    skip_if_not_installed("DOSE")

    result_df <- data.frame(
        ID = c("A", "B", "C"),
        Description = c("DescA", "DescB", "DescC"),
        GeneRatio = c("2/100", "2/100", "2/100"),
        BgRatio = c("10/1000", "10/1000", "10/1000"),
        pvalue = c(0.01, 0.02, 0.03),
        p.adjust = c(0.01, 0.02, 0.03),
        qvalue = c(0.01, 0.02, 0.03),
        geneID = c("x1/x2", "x3/x4", "x5/x6"),
        Count = c(2, 2, 2),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("A", "B", "C")

    termsim <- matrix(0.5, nrow = 3, ncol = 3,
                      dimnames = list(c("DescA", "DescB", "DescC"),
                                      c("DescA", "DescB", "DescC")))
    diag(termsim) <- 1

    mock_result <- new("enrichResult",
                       result = result_df,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,
                       organism = "test",
                       ontology = "test",
                       gene = paste0("x", 1:6),
                       keytype = "SYMBOL",
                       universe = character(0),
                       geneSets = list(
                           A = c("x1", "x2"),
                           B = c("x3", "x4"),
                           C = c("x5", "x6")
                       ),
                       readable = FALSE,
                       termsim = termsim,
                       method = "JC")

    result <- .graphFromEnrichResult(mock_result, c("DescA", "DescC"))

    expect_type(result$geneSet, "list")
    expect_equal(length(result$geneSet), 2)
})

test_that(".graphFromEnrichResult uses custom color parameter", {
    skip_if_not_installed("DOSE")

    result_df <- data.frame(
        ID = c("G1", "G2"),
        Description = c("P1", "P2"),
        GeneRatio = c("2/100", "2/100"),
        BgRatio = c("10/1000", "10/1000"),
        pvalue = c(0.01, 0.02),
        p.adjust = c(0.01, 0.02),
        qvalue = c(0.001, 0.002),
        geneID = c("a1/a2", "a3/a4"),
        Count = c(2, 2),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("G1", "G2")

    termsim <- matrix(c(1, 0.6, 0.6, 1), nrow = 2, ncol = 2,
                      dimnames = list(c("P1", "P2"), c("P1", "P2")))

    mock_result <- new("enrichResult",
                       result = result_df,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,
                       organism = "test",
                       ontology = "test",
                       gene = paste0("a", 1:4),
                       keytype = "SYMBOL",
                       universe = character(0),
                       geneSets = list(
                           G1 = c("a1", "a2"),
                           G2 = c("a3", "a4")
                       ),
                       readable = FALSE,
                       termsim = termsim,
                       method = "JC")

    result <- .graphFromEnrichResult(mock_result, 2, color = "qvalue")

    expect_s3_class(result$graph, "igraph")
    expect_true("color" %in% vertex_attr_names(result$graph))
})

test_that(".graphFromEnrichResult uses custom min_edge parameter", {
    skip_if_not_installed("DOSE")

    result_df <- data.frame(
        ID = c("X1", "X2"),
        Description = c("Proc1", "Proc2"),
        GeneRatio = c("2/100", "2/100"),
        BgRatio = c("10/1000", "10/1000"),
        pvalue = c(0.01, 0.02),
        p.adjust = c(0.01, 0.02),
        qvalue = c(0.01, 0.02),
        geneID = c("m1/m2", "m3/m4"),
        Count = c(2, 2),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("X1", "X2")

    termsim <- matrix(c(1, 0.3, 0.3, 1), nrow = 2, ncol = 2,
                      dimnames = list(c("Proc1", "Proc2"), c("Proc1", "Proc2")))

    mock_result <- new("enrichResult",
                       result = result_df,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,
                       organism = "test",
                       ontology = "test",
                       gene = paste0("m", 1:4),
                       keytype = "SYMBOL",
                       universe = character(0),
                       geneSets = list(
                           X1 = c("m1", "m2"),
                           X2 = c("m3", "m4")
                       ),
                       readable = FALSE,
                       termsim = termsim,
                       method = "JC")

    result <- .graphFromEnrichResult(mock_result, 2, min_edge = 0.5)

    expect_s3_class(result$graph, "igraph")
})

test_that(".graphFromEnrichResult uses custom size_edge parameter", {
    skip_if_not_installed("DOSE")

    result_df <- data.frame(
        ID = c("T1", "T2"),
        Description = c("Term1", "Term2"),
        GeneRatio = c("2/100", "2/100"),
        BgRatio = c("10/1000", "10/1000"),
        pvalue = c(0.01, 0.02),
        p.adjust = c(0.01, 0.02),
        qvalue = c(0.01, 0.02),
        geneID = c("z1/z2", "z3/z4"),
        Count = c(2, 2),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("T1", "T2")

    termsim <- matrix(c(1, 0.7, 0.7, 1), nrow = 2, ncol = 2,
                      dimnames = list(c("Term1", "Term2"), c("Term1", "Term2")))

    mock_result <- new("enrichResult",
                       result = result_df,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,
                       organism = "test",
                       ontology = "test",
                       gene = paste0("z", 1:4),
                       keytype = "SYMBOL",
                       universe = character(0),
                       geneSets = list(
                           T1 = c("z1", "z2"),
                           T2 = c("z3", "z4")
                       ),
                       readable = FALSE,
                       termsim = termsim,
                       method = "JC")

    result <- .graphFromEnrichResult(mock_result, 2, size_edge = 2)

    expect_s3_class(result$graph, "igraph")
})

test_that(".graphFromEnrichResult geneSet contains correct gene sets", {
    skip_if_not_installed("DOSE")

    result_df <- data.frame(
        ID = c("GO:001", "GO:002"),
        Description = c("Set1", "Set2"),
        GeneRatio = c("3/100", "2/100"),
        BgRatio = c("10/1000", "10/1000"),
        pvalue = c(0.01, 0.02),
        p.adjust = c(0.01, 0.02),
        qvalue = c(0.01, 0.02),
        geneID = c("gA/gB/gC", "gD/gE"),
        Count = c(3, 2),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("GO:001", "GO:002")

    termsim <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2,
                      dimnames = list(c("Set1", "Set2"), c("Set1", "Set2")))

    mock_result <- new("enrichResult",
                       result = result_df,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,
                       organism = "test",
                       ontology = "test",
                       gene = c("gA", "gB", "gC", "gD", "gE"),
                       keytype = "SYMBOL",
                       universe = character(0),
                       geneSets = list(
                           "GO:001" = c("gA", "gB", "gC"),
                           "GO:002" = c("gD", "gE")
                       ),
                       readable = FALSE,
                       termsim = termsim,
                       method = "JC")

    result <- .graphFromEnrichResult(mock_result, 2)

    expect_equal(names(result$geneSet), c("Set1", "Set2"))
    expect_equal(result$geneSet$Set1, c("gA", "gB", "gC"))
    expect_equal(result$geneSet$Set2, c("gD", "gE"))
})
