#' @importFrom dplyr filter
#' @importFrom data.table fwrite
#' @importFrom stats sd
#' @importFrom scales rescale
#' @importFrom cowplot ggdraw draw_text
#' @keywords internal

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
                    pathways.all, pathwayLabelStyleDefault, stats, breaklineN){
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
        pathway <- .labelBreak(pn, breaklineN)
        #pathway <- paste0(unlist(stringr::str_split(pn, '_'))[-1], collapse='_')
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

.plotGseaTable <- function(
    pathways, stats, fgseaRes, gseaParam=1, breaklineN=30, fontSize=8,
    colwidths=c(0.8,0.9,0.9,3.6,6), pathwayLabelStyle=NULL, render=NULL,
    headerLabelStyle=NULL, valueStyle=NULL, axisLabelStyle=NULL) {

    pathways.all <- pathways
    pathwayLabelStyleDefault <-
      list(size=fontSize, hjust=0, x=0.05, vjust=.4, lineheight=.8)
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

    pathways <- lapply(pathways, function(p) {
        unname(as.vector(na.omit(match(p, names(statsAdj)))))})

    # fixes #40
    pathways <- pathways[vapply(pathways, length, FUN.VALUE=numeric(1)) > 0]
    NES_all <- fgseaRes %>% dplyr::filter(ID %in% names(pathways)) %>%
        dplyr::pull(NES)
    maxNES <- max(abs(NES_all)) # pmax(abs(NES_all))

    ps <- .get_ps(Pathways=pathways, fgseaRes=fgseaRes, maxNES=maxNES,
            valueStyle=valueStyle, headerLabelStyle=headerLabelStyle,
            pathways.all=pathways.all, breaklineN,
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

.barplot <- function(plot.name, type.sig, topN=10,
                        RES_NES_strings){

    RES_NES_top <- RES_NES_strings %>% dplyr::filter(NES > 0)
    NES_down <- RES_NES_strings %>% dplyr::filter(NES < 0) %>%
        dplyr::arrange(NES)
    RES_NES_top10 <- if(nrow(RES_NES_top) > topN) {
        RES_NES_top %>% dplyr::slice(seq_len(topN))
    } else {
        RES_NES_top
    }
    NES_down10 <- if(nrow(NES_down) > topN) {
        NES_down %>% dplyr::slice(seq_len(topN)) %>% dplyr::arrange(desc(NES))
    } else {
        NES_down
    }
    RES_NES_topdown10 <- rbind(RES_NES_top10, NES_down10) %>%
                                dplyr::select(pathway, NES, eval(type.sig))
    RES_NES_topdown10$Type <- dplyr::if_else(RES_NES_topdown10$NES > 0,
                                            "risk", "protect")
    RES_NES_topdown10$pathway <- factor(RES_NES_topdown10$pathway,
                                        levels=rev(RES_NES_topdown10$pathway))
    tab1 <- RES_NES_topdown10 %>%
                    dplyr::mutate(y=stringr::str_wrap(pathway, width=26))
    p1 <- ggpubr::ggbarplot(data=tab1, y="NES", x="y", color="Type",
                        fill="Type", palette=c("#66C2A5", "#FC8D62"),
                        xlab="", ylab="NES", orientation="horiz",
                        alpha=0.6, sort.val="asc", legend="bottom") +
    ggplot2::geom_vline(xintercept=0, color="#444444") +
    ggplot2::theme(axis.title=ggplot2::element_text(size=7),
                        axis.text.x=ggplot2::element_text(size=7),
                        axis.text.y=ggplot2::element_text(size=5),
                        legend.title=ggplot2::element_text(size=6),
                        legend.text=ggplot2::element_text(size=6))
    #ggplot2::ggsave(filename=plot.name, plot=p1, path=output_path,
    #                width=4, height=6)
    return(p1)
}






