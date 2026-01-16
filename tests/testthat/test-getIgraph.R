library(testthat)
library(igraph)

test_that(".getIgraph extracts graph with numeric nCategory", {
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

    g <- .getIgraph(mock_result, 2, "pvalue", 1, 0.2)

    expect_s3_class(g, "igraph")
    expect_true(vcount(g) <= 2)
})

test_that(".getIgraph extracts graph with character nCategory", {
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

    g <- .getIgraph(mock_result, c("Process 1", "Process 3"), "pvalue", 1, 0.2)

    expect_s3_class(g, "igraph")
})

test_that(".getIgraph errors when nCategory is 0", {
    skip_if_not_installed("DOSE")

    result_df <- data.frame(
        ID = character(0),
        Description = character(0),
        GeneRatio = character(0),
        BgRatio = character(0),
        pvalue = numeric(0),
        p.adjust = numeric(0),
        qvalue = numeric(0),
        geneID = character(0),
        Count = numeric(0),
        stringsAsFactors = FALSE
    )

    termsim <- matrix(numeric(0), nrow = 0, ncol = 0)

    mock_result <- new("enrichResult",
                       result = result_df,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,
                       organism = "test",
                       ontology = "test",
                       gene = character(0),
                       keytype = "SYMBOL",
                       universe = character(0),
                       geneSets = list(),
                       readable = FALSE,
                       termsim = termsim,
                       method = "JC")

    expect_error(
        .getIgraph(mock_result, 0, "pvalue", 1, 0.2),
        "no enriched term found"
    )
})

test_that(".getIgraph passes parameters to buildEmapGraph", {
    skip_if_not_installed("DOSE")

    result_df <- data.frame(
        ID = c("GO:0001", "GO:0002"),
        Description = c("Path 1", "Path 2"),
        GeneRatio = c("2/100", "2/100"),
        BgRatio = c("10/1000", "10/1000"),
        pvalue = c(0.01, 0.02),
        p.adjust = c(0.01, 0.02),
        qvalue = c(0.01, 0.02),
        geneID = c("g1/g2", "g3/g4"),
        Count = c(2, 2),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("GO:0001", "GO:0002")

    termsim <- matrix(c(1.0, 0.6, 0.6, 1.0),
                      nrow = 2, ncol = 2,
                      dimnames = list(c("Path 1", "Path 2"),
                                      c("Path 1", "Path 2")))

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
                           "GO:0001" = c("g1", "g2"),
                           "GO:0002" = c("g3", "g4")
                       ),
                       readable = FALSE,
                       termsim = termsim,
                       method = "JC")

    g <- .getIgraph(mock_result, 2, "qvalue", 2, 0.3)

    expect_s3_class(g, "igraph")
    if (ecount(g) > 0) {
        expect_true("width" %in% edge_attr_names(g))
    }
})

test_that(".getIgraph handles single category", {
    skip_if_not_installed("DOSE")

    result_df <- data.frame(
        ID = "GO:0001",
        Description = "Single Process",
        GeneRatio = "2/100",
        BgRatio = "10/1000",
        pvalue = 0.01,
        p.adjust = 0.01,
        qvalue = 0.01,
        geneID = "gene1/gene2",
        Count = 2,
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- "GO:0001"

    termsim <- matrix(1.0, nrow = 1, ncol = 1,
                      dimnames = list("Single Process", "Single Process"))

    mock_result <- new("enrichResult",
                       result = result_df,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,
                       organism = "test",
                       ontology = "test",
                       gene = c("gene1", "gene2"),
                       keytype = "SYMBOL",
                       universe = character(0),
                       geneSets = list("GO:0001" = c("gene1", "gene2")),
                       readable = FALSE,
                       termsim = termsim,
                       method = "JC")

    g <- .getIgraph(mock_result, 1, "pvalue", 1, 0.2)

    expect_s3_class(g, "igraph")
    expect_equal(vcount(g), 1)
})

test_that(".getIgraph uses termsim and method from object", {
    skip_if_not_installed("DOSE")

    result_df <- data.frame(
        ID = c("ID1", "ID2"),
        Description = c("Desc1", "Desc2"),
        GeneRatio = c("2/100", "2/100"),
        BgRatio = c("10/1000", "10/1000"),
        pvalue = c(0.01, 0.02),
        p.adjust = c(0.01, 0.02),
        qvalue = c(0.01, 0.02),
        geneID = c("g1/g2", "g3/g4"),
        Count = c(2, 2),
        stringsAsFactors = FALSE
    )
    rownames(result_df) <- c("ID1", "ID2")

    termsim <- matrix(c(1.0, 0.7, 0.7, 1.0),
                      nrow = 2, ncol = 2,
                      dimnames = list(c("Desc1", "Desc2"),
                                      c("Desc1", "Desc2")))

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
                           ID1 = c("g1", "g2"),
                           ID2 = c("g3", "g4")
                       ),
                       readable = FALSE,
                       termsim = termsim,
                       method = "Wang")

    g <- .getIgraph(mock_result, 2, "pvalue", 1, 0.2)

    expect_s3_class(g, "igraph")
})
