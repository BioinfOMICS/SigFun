## internal Function for pathway term line break
.labelBreak  <- function(pathwayName, breaklineN){
    pathwayName <- gsub("_", " ", sub("^[^_]+_", "", pathwayName))
    pathwayBreak <- yulab.utils::str_wrap(pathwayName, breaklineN)
    pathwayBreak <- vapply(
        stringr::str_split(pathwayBreak, "\n"),
        function(x) {
            if(length(x) < 3) { paste(x, collapse=" ") } else {
                paste0(paste0(x[seq_len(3)], collapse=" "), '...') }
        }, character(1))
    pathwayBreak <- yulab.utils::str_wrap(pathwayBreak, breaklineN)
    maxBreak <- which(stringr::str_count(pathwayBreak, "\n") >= 3)
    while(length(maxBreak) != 0 ){
        pathwayBreak[maxBreak] <- vapply(
            stringr::str_split(pathwayBreak[maxBreak], "\n"),
            function(x) {
                paste0(paste0(x[seq_len(3)], collapse=" "), '...')
            }, character(1))
        pathwayBreak[maxBreak] <- yulab.utils::str_wrap(
            pathwayBreak[maxBreak], breaklineN)
        maxBreak <- which(stringr::str_count(pathwayBreak, "\n") >= 3)
    }
    return(pathwayBreak)
}

## Internal function for color scaling
.scale_fill <- function(value, neg="#4169E1", zero="white", pos="#FF5151"){
    if(min(value, na.rm=TRUE) >= 0){
        ggplot2::scale_fill_gradient2(
            low=zero, high=pos, limits=c(0, max(value)), oob=scales::squish)
    }else if(max(value, na.rm=TRUE) <= 0){
        ggplot2::scale_fill_gradient2(
            low=neg, high=zero, limits=c(min(value), 0), oob=scales::squish)
    }else{
        ggplot2::scale_fill_gradient2(
            low=neg, mid=zero, high=pos,
            limits=c(min(value), max(value)), oob=scales::squish)
    }
}

## Internal function for color scaling in bar plot
.scale_bar <- function(value){
    if(min(value) >= 0){
        ggplot2::scale_fill_gradient(high="#F8806A", low="#FF5151")
    }else if(max(value) <= 0){
        ggplot2::scale_fill_gradient(high="#deebf7", low="#4169E1")
    }else{
        ggplot2::scale_fill_gradient2(
            low='#4169E1', mid='white', high='#FF5151')
    }
}

## Internal function for chrod plot
.chrod <- function(data, fontSize) {
    pathways <- unique(data$name)
    genes <- unique(data$Gene)
    all_nodes <- union(pathways, genes)
    color_num <- length(all_nodes)
    grid_colors <- viridisLite::viridis(color_num, option = "D")
    names(grid_colors) <- all_nodes
    circlize::circos.clear()
    circlize::circos.par(
        gap.after = c(rep(1, length(genes) - 1), 10,
                      rep(1, length(pathways) - 1), 10)
    )
    unique_sectors <- length(all_nodes)
    max_gap <- 360 / (unique_sectors * 2)
    suggested_gap <- min(5, max_gap)
    circlize::chordDiagram(
        data,
        annotationTrack = "grid",
        annotationTrackHeight = 0.05,
        big.gap = suggested_gap,
        grid.col = grid_colors,
        preAllocateTracks = list(track1 = list(bg.border = NA, track.height = 0.1))
    )
    circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
        xlim <- circlize::get.cell.meta.data("xlim")
        ylim <- circlize::get.cell.meta.data("ylim")
        sector.name <- circlize::get.cell.meta.data("sector.index")
        if (sector.name %in% pathways) {
            circlize::circos.text(
                mean(xlim), ylim[1] + 0.2, sector.name,
                facing = "inside", niceFacing = TRUE,
                adj = c(0.5, 0), cex = fontSize
            )
        } else if (sector.name %in% genes) {
            circlize::circos.text(
                mean(xlim), ylim[1] + 0.2, sector.name,
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = fontSize
            )
        }
    }, bg.border = NA)
    circlize::circos.clear()
}

## Internal function for upset plot
.upsetData <- function(geneSets, res.gsea, breaklineN, type='box'){
    category <- split(geneSets[,1], geneSets[, 2])
    Rawdata <- tibble::tibble(Description=category, gene=names(category),
                              Coef=res.gsea@geneList[names(category)]) |>
        tidyr::unnest(Description) |> dplyr::mutate(present=TRUE) |>
        tidyr::pivot_wider(
            names_from=Description, values_from=present, values_fill=FALSE) |>
        tidyr::gather(key='categoryID', value='value', -gene, -Coef) |>
        dplyr::mutate(name=.labelBreak(categoryID, breaklineN))
    data <- Rawdata |>
        dplyr::select(gene, name, Coef, value) |>
        dplyr::arrange(name) |> tidyr::spread(key='name', value='value')
    pathway_cols <- data |> dplyr::select(where(is.logical)) |> colnames()
    data <- data |> dplyr::rowwise() |>
        dplyr::mutate(
            Path_Combination=paste(
                names(
                    dplyr::across(dplyr::all_of(pathway_cols)))[
                        dplyr::c_across(dplyr::all_of(pathway_cols))],
                collapse="&")) |>
        dplyr::group_by(Path_Combination) |>
        dplyr::mutate(Count=dplyr::n()) |>
        dplyr::ungroup()
    Rawdata <- Rawdata |> dplyr::filter(value == TRUE) |>
        dplyr::group_by(gene) |>
        dplyr::mutate(original_name=paste0(categoryID, collapse=' & '),
                      label_name=paste0(name, collapse=' & ')) |>
        dplyr::select(original_name, label_name, gene, Coef)
    if(type == 'box'){
        data <- data |>
            dplyr::arrange(Path_Combination) |>
            dplyr::mutate(Path_Combination=
                              factor(Path_Combination, levels=unique(Path_Combination)))
    }else if(type == 'bar'){
        data <- data |>
            dplyr::arrange(-Count) |>
            dplyr::mutate(Path_Combination=
                              factor(Path_Combination, levels=unique(Path_Combination)))
        Rawdata <- Rawdata |> dplyr::group_by(original_name) |>
            dplyr::mutate(Count=dplyr::n(),
                          gene=paste0(gene, collapse='/')) |>
            dplyr::select(-Coef) |> dplyr::distinct()
    }
    return(list(data=data, Rawdata=Rawdata))
}



.upsetBarplot <- function(data, fontSize = NULL, fillColor = NULL){
    data |>
        dplyr::select(Path_Combination, Count) |>
        dplyr::distinct() |>
        ggplot2::ggplot(ggplot2::aes(x=Path_Combination, y=Count)) +
        {
            if (!is.null(fillColor)) {
                ggplot2::geom_col(fill = fillColor)
            } else {
                ggplot2::geom_col()
            }
        } +
        .upsetTheme +
        ggplot2::xlab('') +
        ggplot2::ylab('') +
        ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        ggplot2::theme(
            axis.text.y=ggplot2::element_text(color='black', size= fontSize),
            axis.line = ggplot2::element_line(color = "grey40", linewidth = 0.5),
            axis.ticks.x=ggplot2::element_blank()
        )
}

.upsetBoxplot <- function(data, fontSize = NULL, fillColor = NULL) {
    data |>
        ggplot2::ggplot(ggplot2::aes(x = Path_Combination, y = Coef)) +
        {
            if (!is.null(fillColor)) {
                ggplot2::geom_boxplot(fill = fillColor, color = "grey")
            } else {
                ggplot2::geom_boxplot(color = "grey")
            }
        } +
        ggplot2::geom_jitter(width = 0.2, alpha = 0.2) +
        .upsetTheme +
        ggplot2::xlab('') +
        ggplot2::ylab('') +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(color = 'black', size = fontSize),
            axis.line = ggplot2::element_line(color = "grey40", linewidth = 0.5),
            axis.ticks.x = ggplot2::element_blank()
        )
}

.upsetPlot <- function(data, fontSize = NULL, fillColor = NULL){
    pathway_cols <- data |> dplyr::select(where(is.logical)) |> colnames()
    data <- data |>
        dplyr::select(Path_Combination, dplyr::all_of(pathway_cols)) |>
        tidyr::gather(key='Pathway', value='value', -Path_Combination) |>
        dplyr::arrange(Pathway) |>
        dplyr::mutate(Pathway=factor(Pathway,levels=rev(unique(Pathway))))
    data |>
        ggplot2::ggplot(ggplot2::aes(x=Path_Combination, y=Pathway)) +
        ggplot2::geom_point(ggplot2::aes(color=value), size=3) +
        ggplot2::geom_line(data=function(dat) dat[dat$value, ,drop=FALSE],
                           ggplot2::aes(group=Path_Combination), linewidth=1.2, color = fillColor) +
        ggplot2::scale_color_manual(
            values=c(`TRUE`= fillColor, `FALSE`="#E0E0E0")) +
        ggplot2::guides(color="none", fill="none") +
        .upsetTheme + ggplot2::xlab('') + ggplot2::ylab('') +
        ggplot2::theme(
            panel.background=ggplot2::element_blank(),
            axis.text.y=ggplot2::element_text(color='black', size=fontSize),
            axis.ticks=ggplot2::element_blank())
}

.upsetTheme <- ggplot2::theme(
    panel.background=ggplot2::element_rect(fill='white', color=NA),
    panel.border=ggplot2::element_blank(),
    axis.text.x=ggplot2::element_blank(),
    axis.text.y=ggplot2::element_text(colour='black'),
    axis.title.y=ggplot2::element_blank(),
    axis.title.x=ggplot2::element_blank())

.geneOrder <- function(data, distfun, hclustfun){
    hcData <- data |>
        dplyr::select(categoryID, Gene, Coef) |>
        tidyr::spread(key='categoryID', value='Coef') |>
        tibble::column_to_rownames('Gene') |>
        dplyr::mutate(
            dplyr::across(dplyr::everything(), ~tidyr::replace_na(.x, 0)))
    hcData |>
        dist(method=distfun) |>
        hclust(method=hclustfun) |>
        (\(hc) rownames(hcData)[hc$order])()
}


## Internal function for emap and tree plot
.updateN <- function (x, showCategory)
{
    if (!is.numeric(showCategory)) {
        if (inherits(x, "list")) {
            showCategory <- showCategory[showCategory %in% names(x)]
        }
        else {
            showCategory <- intersect(showCategory, x$Description)
        }
        return(showCategory)
    }
    n <- showCategory
    if (inherits(x, "list")) {
        nn <- length(x)
    }
    else {
        nn <- nrow(x)
    }
    if (nn < n) {
        n <- nn
    }
    return(n)
}

## Internal function for emap plot
.graphFromEnrichResult <- function (x, showCategory = 30, color = "p.adjust", min_edge = 0.2,
                                    size_edge = 0.5)
{
    n <- .updateN(x, showCategory)
    y <- as.data.frame(x)
    g <- .getIgraph(x = x, nCategory = n, color = color, cex_line = size_edge,
                    min_edge = min_edge)
    gs <- .extractGeneSets(x, n)
    return(list(graph = g, geneSet = gs))
}


.getIgraph <- function (x, nCategory, color, cex_line, min_edge)
{
    y <- as.data.frame(x)
    geneSets <- DOSE::geneInCategory(x)
    if (is.numeric(nCategory)) {
        y <- y[1:nCategory, ]
    }
    else {
        y <- y[match(nCategory, y$Description), ]
        nCategory <- length(nCategory)
    }
    if (nCategory == 0) {
        stop("no enriched term found...")
    }
    .buildEmapGraph(enrichDf = y, geneSets = geneSets, color = color,
                    cex_line = cex_line, min_edge = min_edge, pair_sim = x@termsim,
                    method = x@method)
}



.buildEmapGraph <- function (enrichDf, geneSets, color, cex_line, min_edge, pair_sim,
                             method)
{
    if (!is.numeric(min_edge) | min_edge < 0 | min_edge > 1) {
        stop("\"min_edge\" should be a number between 0 and 1.")
    }
    if (is.null(dim(enrichDf)) | nrow(enrichDf) == 1) {
        g <- igraph::make_empty_graph(n = 0, directed = FALSE)
        g <- igraph::add_vertices(g, nv = 1)
        igraph::V(g)$name <- as.character(enrichDf$Description)
        igraph::V(g)$color <- "red"
        return(g)
    }
    else {
        w <- pair_sim[as.character(enrichDf$Description), as.character(enrichDf$Description)]
    }
    wd <- reshape2::melt(w)
    wd <- wd[wd[, 1] != wd[, 2], ]
    wd <- wd[!is.na(wd[, 3]), ]
    if (method != "JC") {
        wd[, 1] <- enrichDf[wd[, 1], "Description"]
        wd[, 2] <- enrichDf[wd[, 2], "Description"]
    }
    g <- igraph::graph_from_data_frame(wd[, -3], directed = FALSE)
    igraph::E(g)$width <- sqrt(wd[, 3] * 5) * cex_line
    igraph::E(g)$weight <- wd[, 3]
    # igraph >= 2.0.0: delete.edges() is deprecated
    g <- igraph::delete_edges(g, igraph::E(g)[wd[, 3] < min_edge])
    idx <- unlist(sapply(igraph::V(g)$name, function(x) which(x == enrichDf$Description)))
    cnt <- sapply(geneSets[idx], length)
    igraph::V(g)$size <- cnt
    if (color %in% names(enrichDf)) {
        colVar <- enrichDf[idx, color]
    }
    else {
        colVar <- color
    }
    igraph::V(g)$color <- colVar
    return(g)
}


.extractGeneSets <- function (x, n)
{
    n <- .updateN(x, n)
    if (inherits(x, "list")) {
        geneSets <- x
    }
    else {
        geneSets <- DOSE::geneInCategory(x)
        y <- as.data.frame(x)
        geneSets <- geneSets[y$ID]
        names(geneSets) <- y$Description
    }
    if (is.numeric(n)) {
        return(geneSets[1:n])
    }
    return(geneSets[n])
}

## Internal function for tree plot
.fillTermsim <- function (x, keep)
{
    termsim <- x@termsim[keep, keep]
    termsim[which(is.na(termsim))] <- 0
    termsim2 <- termsim + t(termsim)
    for (i in seq_len(nrow(termsim2))) termsim2[i, i] <- 1
    return(termsim2)
}

.groupTree <- function(hc, clus, d, offsetTiplab, nWords,
                       labelFormatCladelab, labelFormatTiplab, offset,
                       leafFontSize, cladeFontSize,
                       groupColor, extend, hilight,
                       cexCategory, align, dotColor = NULL) {
    group <- count <- NULL
    dat <- data.frame(name = names(clus), cls = paste0("cluster_", as.numeric(clus)))
    grp <- apply(table(dat), 2, function(x) names(x[x == 1]))
    p <- ggtree::ggtree(hc, branch.length = "none")
    noids <- lapply(grp, function(x) {
        unlist(lapply(x, function(i) ggtree::nodeid(p, i)))
    })
    roots <- unlist(lapply(noids, function(x) ggtree::MRCA(p, x)))
    p <- ggtree::groupOTU(p, grp, "group")
    actual_groups <- unique(p$data$group[!is.na(p$data$group)])
    rangeX <- max(p$data$x, na.rm = TRUE) - min(p$data$x, na.rm = TRUE)
    if (inherits(offsetTiplab, "rel")) {
        offsetTiplab <- unclass(offsetTiplab)
        offsetTiplab <- offsetTiplab * 0.5 * max(sqrt(d$count / sum(d$count) * cexCategory))
    }
    if (inherits(offset, "rel")) {
        offset <- unclass(offset)
        offset <- offset * rangeX * 1.2 + offsetTiplab
    }
    pdata <- data.frame(
        name = p$data$label,
        color2 = p$data$group
    )
    pdata <- pdata[!is.na(pdata$name), ]
    clusterColor <- unique(pdata$color2)
    clusterNames <- paste0("cluster_", sort(unique(clus)))
    nColor <- length(clusterNames)
    if (is.null(groupColor)) {
        paletteColors <- viridisLite::viridis(nColor, option = "D")
        names(paletteColors) <- clusterNames
        groupColor <- paletteColors
    } else if (is.null(names(groupColor))) {
        names(groupColor) <- clusterNames[seq_along(groupColor)]
    }
    valid_actual_groups <- actual_groups[!is.na(actual_groups)]
    groupColor_mapped <- setNames(
        groupColor[seq_along(valid_actual_groups)],
        as.character(valid_actual_groups)
    )
    if (any(names(groupColor) %in% valid_actual_groups)) {
        groupColor_mapped <- groupColor[names(groupColor) %in% valid_actual_groups]
    }
    p <- p +
        ggplot2::aes(color = .data$group) +
        ggplot2::scale_color_manual(
            values = groupColor_mapped,
            guide = "none",
            na.value = "grey50",
            drop = FALSE
        )
    p <- ggfun::`%<+%`(p, d)
    if (!is.null(labelFormatTiplab)) {
        labelFuncTiplab <- .defaultLabeller(labelFormatTiplab)
        if (is.function(labelFormatTiplab)) {
            labelFuncTiplab <- labelFormatTiplab
        }
        isTip <- p$data$isTip
        p$data$label[isTip] <- labelFuncTiplab(p$data$label[isTip])
    }
    p <- .addCladelab(
        p = p,
        nWords = nWords,
        label_format_cladelab = labelFormatCladelab,
        offset = offset,
        roots = roots,
        fontsize = cladeFontSize,
        group_color = groupColor,
        cluster_color = clusterColor,
        pdata = pdata,
        extend = extend,
        hilight = hilight,
        align = align
    )
    p <- p + ggnewscale::new_scale_colour()
    if (dotColor == "NES") {
        p <- p +
            ggtree::geom_tippoint(ggplot2::aes(color = .data$color, size = .data$count)) +
            ggplot2::scale_color_gradient2(
                low = "#08306b",
                mid = "white",
                high = "#b3200a",
                midpoint = 0,
                name = "NES",
                labels = function(x) sprintf("%.2f", x)
            )
    } else {
        p <- p +
            ggtree::geom_tippoint(ggplot2::aes(color = -log10(.data$color), size = .data$count)) +
            ggplot2::scale_color_gradient(
                low = "#fcae91",
                high = "#b3200a",
                name = bquote(-log[10](.(dotColor))),
                labels = function(x) sprintf("%.1f", x)
            )
    }
    p <- p + ggtree::geom_tiplab(
        offset = offsetTiplab,
        hjust = 0,
        show.legend = FALSE,
        align = FALSE,
        linewidth = 0,
        size = leafFontSize
    )
    return(p)
}

.defaultLabeller <- function (n)
{
    fun <- function(str) {
        str <- gsub("_", " ", str)
        yulab.utils::str_wrap(str, n)
    }
    structure(fun, class = "labeller")
}

.addCladelab <- function (p, nWords, label_format_cladelab, offset, roots, fontsize,
                          group_color, cluster_color, pdata, extend, hilight, align)
{
    cluster_label <- sapply(cluster_color, .getWordcloud, ggData = pdata,
                            nWords = nWords)
    label_func_cladelab <- .defaultLabeller(label_format_cladelab)
    if (is.function(label_format_cladelab)) {
        label_func_cladelab <- label_format_cladelab
    }
    cluster_label <- label_func_cladelab(cluster_label)
    n_color <- length(levels(cluster_color)) - length(cluster_color)
    if (is.null(group_color)) {
        rlang::check_installed("scales", "for `add_cladelab()`.")
        color2 <- (scales::hue_pal())(length(roots) + n_color)
        if (n_color > 0)
            color2 <- color2[-seq_len(n_color)]
    }
    else {
        color2 <- group_color
    }
    df <- data.frame(node = as.numeric(roots), labels = cluster_label,
                     cluster = cluster_color, color = color2)
    p <- p + ggnewscale::new_scale_colour() + ggtree::geom_cladelab(data = df,
                                                                    mapping = ggplot2::aes(node = node, label = labels, color = cluster),
                                                                    textcolor = "black", extend = extend, show.legend = FALSE,
                                                                    fontsize = fontsize, offset = offset) + ggplot2::scale_color_manual(values = df$color,
                                                                                                                                        guide = "none")
    if (hilight) {
        p <- p + ggtree::geom_hilight(data = df,
                                      mapping = ggplot2::aes(node = node, fill = cluster),
                                      show.legend = FALSE, align = align) +
            ggplot2::scale_fill_manual(values = df$color, guide = "none")
    }
    return(p)
}

.getWordcloud <- function (cluster, ggData, nWords)
{
    `%>%` <- magrittr::`%>%`
    words <- ggData$name %>% gsub(" in ", " ", .) %>% gsub(" [0-9]+ ",
                                                           " ", .) %>% gsub("^[0-9]+ ", "", .) %>% gsub(" [0-9]+$",
                                                                                                        "", .) %>% gsub(" [A-Za-z] ", " ", .) %>% gsub("^[A-Za-z] ",
                                                                                                                                                       "", .) %>% gsub(" [A-Za-z]$", "", .) %>% gsub(" / ",
                                                                                                                                                                                                     " ", .) %>% gsub(" and ", " ", .) %>% gsub(" of ", " ",
                                                                                                                                                                                                                                                .) %>% gsub(",", " ", .) %>% gsub(" - ", " ", .)
    net_tot <- length(words)
    clusters <- unique(ggData$color2)
    words_i <- words[which(ggData$color2 == cluster)]
    sel_tot <- length(words_i)
    sel_w <- .getWordFreq(words_i)
    net_w_all <- .getWordFreq(words)
    net_w <- net_w_all[names(sel_w)]
    tag_size <- (sel_w/sel_tot)/(net_w/net_tot)
    tag_size <- tag_size[order(tag_size, decreasing = TRUE)]
    nWords <- min(nWords, length(tag_size))
    tag <- names(tag_size[seq_len(nWords)])
    dada <- strsplit(words_i, " ")
    len <- vapply(dada, length, FUN.VALUE = 1)
    rank <- NULL
    for (i in seq_len(sel_tot)) {
        rank <- c(rank, seq_len(len[i]))
    }
    word_data <- data.frame(word = unlist(dada), rank = rank)
    word_rank1 <- stats::aggregate(rank ~ word, data = word_data,
                                   sum)
    rownames(word_rank1) <- word_rank1[, 1]
    word_rank1 <- word_rank1[names(sel_w), ]
    word_rank1[, 2] <- word_rank1[, 2]/as.numeric(sel_w)
    tag_order <- word_rank1[tag, ]
    tag_order <- tag_order[order(tag_order[, 2]), ]
    tag_clu_i <- paste(tag_order$word, collapse = " ")
}

.getWordFreq <- function (wordd)
{
    dada <- strsplit(wordd, " ")
    didi <- table(unlist(dada))
    didi <- didi[order(didi, decreasing = TRUE)]
    word_name <- names(didi)
    fun_num_w <- function(ww) {
        sum(vapply(dada, function(w) {
            ww %in% w
        }, FUN.VALUE = 1))
    }
    word_num <- vapply(word_name, fun_num_w, FUN.VALUE = 1)
    word_w <- word_num[order(word_num, decreasing = TRUE)]
}

## Internal function for gsea plot
.getGsdata <- function (x, geneSetID)
{
    if (length(geneSetID) == 1) {
        gsdata <- .gsInfo(x, geneSetID)
        return(gsdata)
    }
    yulab.utils::rbindlist(lapply(geneSetID, .gsInfo, object = x))
}

.gsInfo <- function (object, geneSetID)
{
    geneList <- object@geneList
    if (is.numeric(geneSetID))
        geneSetID <- object@result[geneSetID, "ID"]
    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- .gseaScores(geneList, geneSet, exponent, fortify = TRUE)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore))/20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList
    if (length(object@gene2Symbol) == 0) {
        df$gene <- names(geneList)
    }
    else {
        df$gene <- object@gene2Symbol[names(geneList)]
    }
    df$Description <- object@result[geneSetID, "Description"]
    return(df)
}

.gseaScores <- function (geneList, geneSet, exponent = 1, fortify = FALSE)
{
    geneSet <- intersect(geneSet, names(geneList))
    N <- length(geneList)
    Nh <- length(geneSet)
    Phit <- Pmiss <- numeric(N)
    hits <- names(geneList) %in% geneSet
    Phit[hits] <- abs(geneList[hits])^exponent
    NR <- sum(Phit)
    Phit <- cumsum(Phit/NR)
    Pmiss[!hits] <- 1/(N - Nh)
    Pmiss <- cumsum(Pmiss)
    runningES <- Phit - Pmiss
    max.ES <- max(runningES)
    min.ES <- min(runningES)
    if (abs(max.ES) > abs(min.ES)) {
        ES <- max.ES
    }
    else {
        ES <- min.ES
    }
    df <- data.frame(x = seq_along(runningES), runningScore = runningES,
                     position = as.integer(hits))
    if (fortify == TRUE) {
        return(df)
    }
    df$gene = names(geneList)
    res <- list(ES = ES, runningES = df)
    return(res)
}

.tableGrob2 <- function (d, p = NULL)
{
    d <- d[order(rownames(d)), ]
    rlang::check_installed("gridExtra", "for `tableGrob2()`.")
    tp <- gridExtra::tableGrob(d)
    if (is.null(p)) {
        return(tp)
    }
    p_data <- ggplot_build(p)$data[[1]]
    p_data <- p_data[order(p_data[["group"]]), ]
    pcol <- unique(p_data[["colour"]])
    j <- which(tp$layout$name == "rowhead-fg")
    for (i in seq_along(pcol)) {
        tp$grobs[j][[i + 1]][["gp"]] <- gpar(col = pcol[i])
    }
    return(tp)
}

## Internal function for cnet plot
.subsetCnetList <- function (x, showCategory)
{
    if (!is.numeric(showCategory)) {
        return(x[names(x) %in% showCategory])
    }
    n <- length(x)
    if (length(showCategory) == 1) {
        showCategory <- seq(showCategory)
    }
    if (any(showCategory) > n) {
        msg <- sprintf("any showCategory value that is large than %d will be removed.",
                       n)
        message(msg)
    }
    showCategory <- showCategory[showCategory <= n]
    return(x[showCategory])
}

.subsetCnetListItem <- function (x, showItem = "all")
{
    if (length(showItem) == 1 && showItem == "all")
        return(x)
    lapply(x, function(y) y[y %in% showItem])
}


.list2graph <- function (inputList, directed = FALSE)
{
    x <- .list2df(inputList)
    g <- igraph::graph_from_data_frame(x, directed = directed)
    igraph::V(g)$.isCategory <- igraph::V(g)$name %in% names(inputList)
    size <- vapply(inputList, length, FUN.VALUE = numeric(1))
    igraph::V(g)$size <- igraph::degree(g)
    return(g)
}


.getEdgeData <- function (g)
{
    e <- as.data.frame(igraph::as_edgelist(g))
    enames <- igraph::edge_attr_names(g)
    if (length(enames) > 0) {
        for (eattr in enames) {
            e[[eattr]] <- igraph::edge_attr(g, eattr)
        }
    }
    return(e)
}



