library(testthat)

testthat::test_that(".addCladelab returns modified plot object", {
    skip_if_not_installed("ggtree")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggnewscale")
    skip_if_not_installed("ape")

    tree <- ape::read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
    p <- ggtree::ggtree(tree)

    pdata <- data.frame(
        name = c("gene1 pathway", "gene2 pathway", "gene3 system", "gene4 system"),
        color2 = c("cluster1", "cluster1", "cluster2", "cluster2"),
        stringsAsFactors = FALSE
    )

    cluster_color <- factor(c("cluster1", "cluster2"))

    result <- .addCladelab(
        p = p,
        nWords = 2,
        label_format_cladelab = 30,
        offset = 0.5,
        roots = c(5, 6),
        fontsize = 3,
        group_color = c("red", "blue"),
        cluster_color = cluster_color,
        pdata = pdata,
        extend = 0.5,
        hilight = FALSE,
        align = "left"
    )

    testthat::expect_s3_class(result, "gg")
    testthat::expect_s3_class(result, "ggplot")
})

testthat::test_that(".addCladelab works with hilight = TRUE", {
    skip_if_not_installed("ggtree")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ape")

    tree <- ape::read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
    p <- ggtree::ggtree(tree)

    pdata <- data.frame(
        name = c("term1 alpha", "term2 beta", "term3 gamma", "term4 delta"),
        color2 = c("X", "X", "Y", "Y"),
        stringsAsFactors = FALSE
    )

    result <- .addCladelab(
        p = p,
        nWords = 1,
        label_format_cladelab = 15,
        offset = 0.2,
        roots = c(5, 6),
        fontsize = 3,
        group_color = c("green", "orange"),
        cluster_color = factor(c("X", "Y")),
        pdata = pdata,
        extend = 0.4,
        hilight = TRUE,
        align = "left"
    )

    testthat::expect_s3_class(result, "ggplot")
})

testthat::test_that(".addCladelab works with custom label function", {
    skip_if_not_installed("ggtree")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ape")

    tree <- ape::read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
    p <- ggtree::ggtree(tree)

    pdata <- data.frame(
        name = c("apple banana", "cherry date", "elderberry fig", "grape honeydew"),
        color2 = c("G1", "G1", "G2", "G2"),
        stringsAsFactors = FALSE
    )

    custom_labeller <- function(x) toupper(x)

    result <- .addCladelab(
        p = p,
        nWords = 1,
        label_format_cladelab = custom_labeller,
        offset = 0.5,
        roots = c(5, 6),
        fontsize = 3,
        group_color = c("purple", "cyan"),
        cluster_color = factor(c("G1", "G2")),
        pdata = pdata,
        extend = 0.5,
        hilight = FALSE,
        align = "left"
    )

    testthat::expect_s3_class(result, "ggplot")
})
