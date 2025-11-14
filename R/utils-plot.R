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
    p <- p %<+% d
    if (!is.null(labelFormatTiplab)) {
        labelFuncTiplab <- enrichplot:::default_labeller(labelFormatTiplab)
        if (is.function(labelFormatTiplab)) {
            labelFuncTiplab <- labelFormatTiplab
        }
        isTip <- p$data$isTip
        p$data$label[isTip] <- labelFuncTiplab(p$data$label[isTip])
    }
    p <- enrichplot:::add_cladelab(
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
