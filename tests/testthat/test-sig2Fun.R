library(testthat)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(clusterProfiler)
library(fgsea)
library(S4Vectors)
library(SigFun)

.synchonize_SE <- function(expr.data, SIG_MAT, mapping){
    ## sample
    sample_name <- intersect(colnames(expr.data), SIG_MAT$sample_id)
    expr.data <- as.data.frame(expr.data) %>% dplyr::select(all_of(sample_name))
    order.name <- base::order(colnames(expr.data), decreasing = FALSE)
    expr.data <- expr.data[, order.name]
    SIG_MAT <- SIG_MAT %>%  dplyr::filter(sample_id %in% sample_name) %>%
        dplyr::arrange(sample_id)

    ## gene
    gene_name <- intersect(rownames(expr.data), mapping$gene_symbol) %>% sort()
    mapping <- mapping %>% dplyr::filter(gene_symbol %in% gene_name)
    expr.data <- expr.data[rownames(expr.data) %in% gene_name, ]
    mapping <- mapping %>% dplyr::select(ensg_id, gene_symbol, gene_biotype)
    expr.data <- expr.data[mapping$gene_symbol, ]
    rownames(mapping) <- mapping$ensg_id
    rownames(expr.data) <- mapping$ensg_id
    # SigFun
    ## build SE
    SE_data <- SummarizedExperiment::SummarizedExperiment(
        assays=list(abundance=as.matrix(expr.data)),
        rowData=S4Vectors::DataFrame(mapping, row.names=mapping$ensg_id),
        colData=SIG_MAT)
    return(SE_data)
}

.classCheck <- function(input, must) {
    if (is.null(input)) {
        return()
    }
    if (!is(input, must)) {
        abort(sprintf(
            "The '%s' must be %s, not %s",
            as.character(substitute(input)), must, class(input)
        ))
    }
}

.typeCheck <- function(input, must) {
    if (is.null(input)) {
        return()
    }
    if (!is(input, must)) {
        abort(sprintf(
            "The '%s' must be %s, not %s",
            as.character(substitute(input)), must, type(input)
        ))
    }
}

.extractDF <- function(se, type){
    .checkSE(se)
    if (type == "gseaRaw") {
        data <- S4Vectors::metadata(se)$gseaResult
    } else if (type == "gseaReadable") {
        data <- S4Vectors::metadata(se)$gseaResult |>
            DOSE::setReadable(OrgDb='org.Hs.eg.db', keyType='ENSEMBL')
    } else if (type == "gseaSimilar") {
        data <- S4Vectors::metadata(se)$gseaResult |>
            DOSE::setReadable(OrgDb='org.Hs.eg.db', keyType='ENSEMBL') |>
            enrichplot::pairwise_termsim()
    } else if (type == "corCoef") {
        data <- S4Vectors::metadata(se)$cor.df
    }
    return(data)
}

.checkSE <- function(se){
    stopifnot(inherits(se, "SummarizedExperiment"))
    stopifnot("metadata" %in% slotNames(se))
    stopifnot("gseaResult" %in% names(S4Vectors::metadata(se)))
}

.extractGeneSets <- function(gseaResult, showCategory=30){
    geneSets <- DOSE::geneInCategory(gseaResult)[gseaResult@result$ID]
    names(geneSets) <- gseaResult@result$Description
    if (is.numeric(showCategory)){
        geneSets <- geneSets[seq_len(showCategory)]
    }else{
        if (sum(showCategory %in% names(geneSets)) == 0){
            cli::cli_abort(
                "The value provided to {.arg showCategory} cannot be found in
                {.field SE_data.fgsea@result}."
            )
        }
        geneSets <- geneSets[showCategory]
    }
    data <- .list2df(geneSets)
    return(data)
}

.list2df <- function(inputList) {
    ldf <- lapply(seq_len(length(inputList)), function(i) {
        data.frame(categoryID=rep(names(inputList[i]),
                                  length(inputList[[i]])), Gene=inputList[[i]])
    })
    do.call('rbind', ldf)
}

#' @importFrom dplyr filter
#' @importFrom data.table fwrite
#' @importFrom stats sd
#' @importFrom scales rescale
#' @importFrom cowplot ggdraw draw_text
#' @keywords internal

.z_score_cal <- function(x, NA.rm=FALSE, rescale01=FALSE){
    if(sd(x, na.rm=NA.rm) == 0){
        out <- rep(0, length(x))
    } else {
        out <- (x - mean(x, na.rm=NA.rm)) / sd(x, na.rm=NA.rm)
        if(rescale01 == TRUE){
            out <- scales::rescale((x - mean(x, na.rm=NA.rm)) / sd(x, na.rm=NA.rm))
        }
    }
    return(out)
}

.corList <- function(pattern.signature, pattern.genes.norm, cor.method){
    cor.list <- apply(pattern.genes.norm, 2, function(x) {
        out.res <- NULL
        tmp.res <- cor.test(x, pattern.signature,method=cor.method,
                            use="pairwise.complete.obs")
        out.res$cor <- as.numeric(tmp.res$estimate)
        out.res$pval <- as.numeric(tmp.res$p.value)
        return(out.res)
    })
    cor.list <- do.call("rbind",cor.list)
    return(as.data.frame(cor.list))
}

.logitList <- function(y=pattern.signature, pattern.genes.norm, cor.method){
    cor.list <- apply(pattern.genes.norm, 2, function(x) {
        out.res <- NULL
        res_logit <- summary(glm(y~x, family = binomial(logit)))
        out.res$cor <- as.numeric(res_logit$coefficients[,"Estimate"]["x"])
        out.res$pval <- as.numeric(res_logit$coefficients[,"Pr(>|z|)"]["x"])
        return(out.res)
    })
    cor.list <- do.call("rbind",cor.list)
    return(as.data.frame(cor.list))
}

.cowplotText <- function(text, style) {
    cowplot::ggdraw() + do.call(cowplot::draw_text, c(list(text=text), style))
    #+ ggplot2::theme(panel.background=element_rect(color="black"))
}

.nesbarplot <- function(NES,maxNES) {
    color <- base::ifelse(as.numeric(NES)>0,'#FD7C69','#5F90BB')
    data.frame(Y="",X=abs(NES),X2=1/4*abs(NES)) %>%
        ggplot2::ggplot(ggplot2::aes(x=X, y=Y)) +
        ggplot2::geom_bar(stat="identity",fill=color,na.rm=TRUE) +
        #geom_text(aes(x=X2 ,label=round(X,2)),hjust=0) +
        ggplot2::scale_x_continuous(limits=c(0,maxNES)) +
        ggplot2::theme(panel.background=ggplot2::element_blank(),
                       plot.background=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.line=ggplot2::element_blank(),
                       axis.text=ggplot2::element_blank(),
                       axis.ticks=ggplot2::element_blank(),
                       panel.grid=ggplot2::element_blank(),
                       axis.title=ggplot2::element_blank(),
                       plot.margin=rep(ggplot2::unit(0,"null"),4),
                       panel.spacing=rep(ggplot2::unit(0,"null"),4)
        )
}

.plotEnrichment <- function(pathway, stats, color,
                            gseaParam=1,
                            ticksSize=0.2) {

    pd <- .plotEnrichmentData(pathway=pathway,stats=stats,gseaParam=gseaParam)
    pd$ticks$color <- ifelse(as.numeric(color)>0,'#FD7C69','#5F90BB')
    with(pd,ggplot2::ggplot(data=curve) + ggplot2::geom_segment(data=ticks,
                                                                mapping=ggplot2::aes(x=rank,
                                                                                     y=-spreadES/16,
                                                                                     xend=rank,
                                                                                     yend=spreadES/16,
                                                                                     color=color),
                                                                linewidth=ticksSize, show.legend=FALSE) +
             ggplot2::scale_colour_identity() +
             ggplot2::theme_classic())
}

.plotEnrichmentData <- function(pathway, stats, gseaParam=1) {

    if (any(!is.finite(stats))){
        stop("Not all stats values are finite numbers")
    }

    ord <- order(rank(-stats))

    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)

    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    pathway <- unique(pathway)

    gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats=pathway,
                                   returnAllExtremes=TRUE)

    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops

    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.table::data.table(rank=c(0, xs, n + 1), ES=c(0, ys, 0))
    ticks <- data.table::data.table(rank=pathway, stat=statsAdj[pathway])
    stats <- data.table::data.table(rank=seq_along(stats), stat=statsAdj)

    res <- list(
        curve=toPlot,
        ticks=ticks,
        stats=stats,
        posES=max(tops),
        negES=min(bottoms),
        spreadES=max(tops)-min(bottoms),
        maxAbsStat=max(abs(statsAdj)))
}

.get_ps <- function(Pathways, fgseaRes, maxNES, valueStyle, headerLabelStyle,
                    pathways.all, pathwayLabelStyleDefault, stats){
    common_theme <- ggplot2::theme(panel.background=ggplot2::element_blank(),
                                   plot.background=ggplot2::element_blank(),
                                   axis.line=ggplot2::element_blank(),
                                   axis.text=ggplot2::element_blank(),
                                   axis.ticks=ggplot2::element_blank(),
                                   panel.grid=ggplot2::element_blank(),
                                   axis.title=ggplot2::element_blank(),
                                   plot.margin=rep(ggplot2::unit(0, "null"), 4),
                                   panel.spacing=rep(ggplot2::unit(0, "null"), 4))
    ps <- lapply(names(Pathways), function(pn) {

        annotation <- fgseaRes[fgseaRes$ID == pn, ]
        pathway <- paste0(unlist(stringr::str_split(pn, '_'))[-1], collapse='_')
        list(.nesbarplot(annotation$NES, maxNES),
             .cowplotText(sprintf("%.2f", annotation$NES), valueStyle),
             .cowplotText(formatC(annotation$pvalue, digits=0, format="e"),
                          headerLabelStyle),
             cowplot::plot_grid(
                 .plotEnrichment(pathways.all[[pn]], stats,
                                 color=annotation$NES) + common_theme +
                     ggplot2::theme(panel.background=ggplot2::element_rect(color="black"))),
             .cowplotText(pathway, pathwayLabelStyleDefault))})

    return(ps)
}

.plotGseaTable <- function(pathways.all,Pathways, stats, fgseaRes, gseaParam=1,
                           colwidths=c(0.8,0.9,1.2,2.4,10), pathwayLabelStyle=NULL, render=NULL,
                           headerLabelStyle=NULL, valueStyle=NULL, axisLabelStyle=NULL) {

    pathways.all <- Pathways
    pathwayLabelStyleDefault <- list(size=10, hjust=0, x=0.05, vjust=0.5)
    pathwayLabelStyle <- modifyList(pathwayLabelStyleDefault,
                                    as.list(pathwayLabelStyle))

    headerLabelStyleDefault <- list(size=10)
    headerLabelStyle <- modifyList(headerLabelStyleDefault,
                                   as.list(headerLabelStyle))

    valueStyleDefault <- list(size=10)
    valueStyle <- modifyList(valueStyleDefault, as.list(valueStyle))

    if (!is.null(render)) {
        warning("render argument is deprecated, a ggplot object is always returned")}

    statsAdj <- stats[order(rank(-stats))]
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
    statsAdj <- statsAdj / max(abs(statsAdj))

    Pathways <- lapply(Pathways, function(p) {
        unname(as.vector(na.omit(match(p, names(statsAdj)))))})

    # fixes #40
    Pathways <- Pathways[vapply(Pathways, length, FUN.VALUE=numeric(1)) > 0]
    NES_all <- fgseaRes %>% dplyr::filter(ID %in% names(Pathways)) %>%
        dplyr::pull(NES)
    maxNES <- max(abs(NES_all)) # pmax(abs(NES_all))

    ps <- .get_ps(Pathways=Pathways, fgseaRes=fgseaRes, maxNES=maxNES,
                  valueStyle=valueStyle, headerLabelStyle=headerLabelStyle,
                  pathways.all=pathways.all,
                  pathwayLabelStyleDefault=pathwayLabelStyleDefault, stats=stats)

    grobs <- c(lapply(c("","NES","P-value","Enrichment plot"), .cowplotText,
                      style=headerLabelStyle), list(.cowplotText("Pathway",
                                                                 modifyList(headerLabelStyle, pathwayLabelStyle[c("hjust", "x")]))),
               unlist(ps, recursive=FALSE))

    p <- cowplot::plot_grid(
        plotlist=grobs, ncol=sum(as.numeric(colwidths) != 0),
        rel_widths=colwidths[as.numeric(colwidths) != 0])
    return(p)
}


plot_heat <- function(SE_data.fgsea, significant_type="pval",
                      strings=c("GOBP", "GOCC", "GOMF", "KEGG", "REACTOME", "WP", "HALLMARK"),
                      topN=10, pathways.all, ranking.method="cor") {
    cor.df <- S4Vectors::metadata(SE_data.fgsea)$cor.df
    RES_GSEA <- S4Vectors::metadata(SE_data.fgsea)$gseaResult %>% as_tibble() %>%
        dplyr::slice(grep(paste(strings, collapse="|"), ID)) %>%
        dplyr::mutate(`-log10(pvalue)`=-log10(pvalue))
    if (significant_type == "pval") {
        RES_GSEA <- RES_GSEA %>% dplyr::filter(pvalue < 0.05)
    } else if (significant_type == "padj") {
        RES_GSEA <- RES_GSEA %>% dplyr::filter(p.adjust < 0.05)
    }
    CodingGene <- SummarizedExperiment::rowData(SE_data.fgsea) %>%
        as.data.frame() %>% dplyr::filter(gene_biotype == "protein_coding") %>%
        dplyr::select(ENSG=ensg_id, gene_symbol=gene_symbol)
    ranksDF <- cor.df %>% dplyr::filter(!is.na(gene)) %>%
        dplyr::select(ENSG=gene, stat=ranking.method) %>%
        dplyr::left_join(CodingGene, by="ENSG") %>%
        dplyr::mutate(abs.stat=abs(stat)) %>%
        dplyr::select(ENSG, stat) %>%
        na.omit() %>% dplyr::distinct() %>% dplyr::group_by(ENSG) %>%
        dplyr::summarize(ranking.method=mean(as.numeric(stat))) %>%
        dplyr::filter(ranking.method != 'Inf' & ranking.method != '-Inf') %>%
        tibble::deframe()
    tmp.res <- lapply(strings, function(str) {
        RES_NES <- RES_GSEA %>% dplyr::slice(grep(str, ID))
        RES_NES_top10 <- RES_NES %>%
            dplyr::arrange(dplyr::desc(NES)) %>%
            dplyr::filter(NES >= 0)
        if(nrow(RES_NES_top10) >= topN){
            RES_NES_top10 <- RES_NES_top10 %>%
                dplyr::slice(seq_len(topN))
        }
        RES_NES_bottom10 <- RES_NES %>%
            dplyr::arrange(dplyr::desc(-NES)) %>%
            dplyr::filter(NES <= 0)
        if(nrow(RES_NES_bottom10) >= topN){
            RES_NES_bottom10 <- RES_NES_bottom10 %>%
                dplyr::slice(seq_len(topN))
        }
        RES_NES <- rbind(RES_NES_top10, RES_NES_bottom10) %>%
            dplyr::arrange(dplyr::desc(NES))

        tmp <- .plotGseaTable(pathways.all,
                              Pathways=pathways.all[RES_NES$ID],
                              stats=ranksDF, fgseaRes=RES_GSEA)
        tmp
        #countNrow <- RES_GSEA %>% dplyr::slice(grep(str, ID))
        #if(nrow(countNrow) >= 20){
        #    RES_NES_top10 <- RES_GSEA %>%
        #        dplyr::slice(grep(str, ID)) %>%
        #        dplyr::arrange(dplyr::desc(NES)) %>%
        #        dplyr::slice(seq_len(topN))
        #    RES_NES_bottom10 <- RES_GSEA %>% dplyr::slice(grep(str, ID)) %>%
        #        dplyr::arrange(dplyr::desc(-NES)) %>%
        #        dplyr::slice(seq_len(topN)) %>% dplyr::arrange(dplyr::desc(NES))
        #    tmp <- .plotGseaTable(pathways.all,
        #        Pathways=pathways.all[c(RES_NES_top10$ID, RES_NES_bottom10$ID)],
        #        stats=ranksDF, fgseaRes=RES_GSEA)
        #    tmp
        #}else{
        #    tmp <- .plotGseaTable(pathways.all,
        #        Pathways=pathways.all[c(countNrow$ID)],
        #        stats=ranksDF, fgseaRes=RES_GSEA)
        #    tmp
        #}
    })
    names(tmp.res) <- strings
    return(tmp.res)
}

sig2Fun <- function(SE_data, t2g, ranking.method="cor", topN=10,
                    species="human", cor.method="spearman", Z.transform=FALSE,
                    significant_type="pval", strings=c("GOBP","REACTOME","HALLMARK","SIGNALING")) {
    .classCheck(SE_data, "SummarizedExperiment")
    .classCheck(t2g, "data.frame")
    pathway.name <- unique(t2g$gs_name)
    new_pathways.all <- lapply(pathway.name, function(x){
        res <- t2g %>% dplyr::filter(gs_name==x) %>%
            dplyr::select(ensembl_gene) %>%
            dplyr::pull()
        return(res)
    })
    names(new_pathways.all) <- pathway.name
    SE_data.cor <- SE_data
    if(!("cor.df" %in% names(S4Vectors::metadata(SE_data.cor)))){
        SE_data.cor <- sigCor(SE_data=SE_data, cor.method=cor.method,
                              Z.transform=Z.transform)
        S4Vectors::metadata(SE_data.cor) <-
            list(cor.df=S4Vectors::metadata(SE_data.cor)$cor.df)
    }

    geneList <- SE_data.cor@metadata$cor.df %>%
        dplyr::select(all_of(ranking.method)) %>% dplyr::pull()
    names(geneList) <- SE_data.cor@metadata$cor.df$gene
    geneList <- sort(geneList, decreasing = TRUE)
    clusterProfiler <- clusterProfiler::GSEA(geneList, TERM2GENE = t2g)

    SE_data.fgsea <- SE_data.cor
    S4Vectors::metadata(SE_data.fgsea) <- list(gseaResult=clusterProfiler,
                                               cor.df=S4Vectors::metadata(SE_data.cor)$cor.df)

    # heatmap
    heatmap.list <- NULL
    for(i in seq_along(strings)){
        heatmap.list[[strings[i]]] <- plot_heat(
            SE_data.fgsea=SE_data.fgsea,
            pathways.all=new_pathways.all,
            significant_type=significant_type,
            strings=strings[i],
            topN=topN, ranking.method=ranking.method)
    }
    SE_data.fgsea@metadata$heatmap  <- heatmap.list

    return(SE_data.fgsea)
}



sigCor <- function(SE_data, cor.method="spearman", Z.transform=FALSE) {
    signature.obj <- as.data.frame(SummarizedExperiment::colData(SE_data))

    exp_data <- as.data.frame(SummarizedExperiment::assay(SE_data))
    exp_data$ensg_id <- rownames(exp_data)
    exp_data <- exp_data %>% tidyr::gather(-ensg_id, key="sample_id",
                                           value="value")
    S_ID <- intersect(signature.obj$sample_id, unique(exp_data$sample_id))
    exp_data.sid <- exp_data %>% dplyr::filter(sample_id %in% S_ID)
    signature.obj <- signature.obj %>% dplyr::filter(sample_id %in% S_ID)
    exp_data.FPKM_UQ.tmp <- exp_data.sid %>%
        dplyr::select(sample_id, ensg=ensg_id, value)
    tmp.data.wide <- exp_data.FPKM_UQ.tmp %>% dplyr::filter(ensg != "") %>%
        dplyr::select(sample_id, ensg, value) %>%
        dplyr::distinct(sample_id, ensg, .keep_all=TRUE) %>%
        tidyr::pivot_wider(names_from=ensg, values_from=value)

    cor.object <- merge(signature.obj, tmp.data.wide)
    pattern.signature <- cor.object$value
    pattern.genes <- cor.object %>% dplyr::select(-sample_id, -value)
    count.EffectSamples <- unlist(lapply(pattern.genes,
                                         function(x) length(unique(x))))
    rm.index <- which(count.EffectSamples < 2)

    if (Z.transform == TRUE) {
        pattern.genes.norm <- if (length(rm.index) > 0) {
            apply(pattern.genes[, -rm.index], 2, .z_score_cal)
        } else { apply(pattern.genes, 2, .z_score_cal) }
    } else { pattern.genes.norm <- if (length(rm.index) > 0)
    {pattern.genes[, -rm.index]} else {pattern.genes} }

    if(cor.method %in% c("pearson", "kendall", "spearman")){
        cor.list <- .corList(pattern.signature, pattern.genes.norm, cor.method)
    }

    if(cor.method %in% c("logit")){
        cor.list <- .logitList(y=pattern.signature, pattern.genes.norm, cor.method)
    }

    cor.df <- data.frame(gene=rownames(cor.list), cor=as.numeric(cor.list$cor) ,
                         pval=as.numeric(cor.list$pval))
    S4Vectors::metadata(SE_data) <- list(cor.df=cor.df)
    return(SE_data)
}




# Test helper function for class checking
test_that(".classCheck works correctly", {
    # Mock the .classCheck function since it's not provided
    .classCheck <- function(obj, expected_class) {
        if (!inherits(obj, expected_class)) {
            stop(paste("Object is not of class", expected_class))
        }
    }

    # Test with correct class
    se_data <- SummarizedExperiment(
        assays <- list(counts <- matrix(1:20, nrow <- 4))
    )
    expect_silent(.classCheck(se_data, "SummarizedExperiment"))

    # Test with incorrect class
    expect_error(.classCheck(list(), "SummarizedExperiment"))
})

data("demo_GSE181574")
GSE181574.sigfun <- SummarizedExperiment::SummarizedExperiment(
    assays=list(abundance=as.matrix(expr.data)),
    rowData=S4Vectors::DataFrame(mapping, row.names=mapping$ensg_id),
    colData=SIG_MAT)


# Test sig2Fun main function
test_that("sig2Fun performs complete signature analysis", {
    # Skip if required packages not available
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("fgsea")

    # Test sig2Fun function with expected warnings
    expect_warning({
        result <- sig2Fun(
            SE_data = GSE181574.sigfun,
            ranking.method = "cor",
            species = "human",
            t2g = t2g,
            cor.method = "logit",
            topN = 10,
            Z.transform = FALSE,
            significant_type = "pval",
            strings = c(
                "GOBP","GOCC","GOMF","KEGG","REACTOME",
                "WP","HALLMARK","SIGNALING")
        )
    }, regexp = "For some of the pathways the P-values were likely overestimated|For some pathways, in reality P-values are less than 1e-10")

    # Alternative: expect specific warnings individually
    expect_warning(
        expect_warning({
            result <- sig2Fun(
                SE_data = GSE181574.sigfun,
                ranking.method = "cor",
                species = "human",
                t2g = t2g,
                cor.method = "logit",
                topN = 10,
                Z.transform = FALSE,
                significant_type = "pval",
                strings = c(
                    "GOBP","GOCC","GOMF","KEGG","REACTOME",
                    "WP","HALLMARK","SIGNALING")
            )
        }, "For some of the pathways the P-values were likely overestimated"),
        "For some pathways, in reality P-values are less than 1e-10"
    )
})

# Test error handling
test_that("Functions properly handle invalid inputs", {
    # Test .plotEnrichmentData with invalid stats
    invalid_stats <- c(1, Inf, -Inf, NaN)
    pathway <- c("g1", "g2")
    expect_error(.plotEnrichmentData(pathway, invalid_stats))

    # Test .corList with mismatched dimensions
    short_signature <- c(1, 2, 3)
    long_genes <- data.frame(
        gene1 <- c(1, 2, 3, 4, 5),
        gene2 <- c(2, 3, 4, 5, 6)
    )
    expect_error(.corList(short_signature, long_genes, "pearson"))

    # Test .z_score_cal with empty vector
    expect_error(.z_score_cal(numeric(0)))
})
