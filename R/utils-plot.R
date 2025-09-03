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
.chrod <- function(data, fontSize){
    color_num <- length(union(unique(data$name), unique(data$Gene)))
    pathways <- unique(data$name)
    genes <- unique(data$Gene)
    circlize::circos.clear()
    circlize::circos.par(
        gap.after=c(rep(1,length(unique(data$Gene))-1),10,
                    rep(1,length(unique(data$name))-1),10))
    unique_sectors <- length(unique(c(data[, 1], data[, 2])))
    max_gap <- 360 / (unique_sectors * 2)
    suggested_gap <- min(5, max_gap)
    p <- circlize::chordDiagram(
        data,
    annotationTrack="grid",  annotationTrackHeight=0.05, big.gap=suggested_gap,
        grid.col=randomcoloR::distinctColorPalette(color_num),
        preAllocateTracks=list(track1=list(bg.border=NA, track.height=0.1)))
    circlize::circos.track(track.index=1, panel.fun=function(x, y) {
        xlim <- circlize::get.cell.meta.data("xlim")
        ylim <- circlize::get.cell.meta.data("ylim")
        sector.name <- circlize::get.cell.meta.data("sector.index")
        if (sector.name %in% pathways) {
            circlize::circos.text(
                mean(xlim), ylim[1]+0.2, sector.name, facing="inside",
                niceFacing=TRUE, adj=c(0.5, 0), cex=fontSize)
        } else if (sector.name %in% genes) {
            circlize::circos.text(
                mean(xlim), ylim[1]+0.2, sector.name, facing="clockwise",
                niceFacing=TRUE, adj=c(0, 0.5), cex=fontSize)
        }
    }, bg.border=NA)
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

.upsetBarplot <- function(data){
    data |> dplyr::select(Path_Combination, Count) |> dplyr::distinct() |>
        ggplot2::ggplot(ggplot2::aes(x=Path_Combination, y=Count)) +
        ggplot2::geom_col() + .upsetTheme + ggplot2::xlab('') +
        ggplot2::ylab('') +
        ggplot2::theme(
            axis.text.y=ggplot2::element_text(color='black', size=12))
}

.upsetBoxplot <- function(data){
    data |>
        ggplot2::ggplot(ggplot2::aes(x=Path_Combination, y=Coef)) +
        ggplot2::geom_boxplot() + ggplot2::geom_jitter(width=.2, alpha=.2) +
        .upsetTheme + ggplot2::xlab('') + ggplot2::ylab('') +
        ggplot2::theme(
            axis.text.y=ggplot2::element_text(color='black', size=12))
}

.upsetPlot <- function(data){
    pathway_cols <- data |> dplyr::select(where(is.logical)) |> colnames()
    data <- data |>
        dplyr::select(Path_Combination, dplyr::all_of(pathway_cols)) |>
        tidyr::gather(key='Pathway', value='value', -Path_Combination) |>
        dplyr::arrange(Pathway) |>
        dplyr::mutate(Pathway=factor(Pathway,levels=rev(unique(Pathway))))
    fontSize <- ifelse(length(unique(data$Pathway)) < 10, 12, 8)
    data |>
        ggplot2::ggplot(ggplot2::aes(x=Path_Combination, y=Pathway)) +
        ggplot2::geom_point(ggplot2::aes(color=value), size=3) +
        ggplot2::geom_line(data=function(dat) dat[dat$value, ,drop=FALSE],
            ggplot2::aes(group=Path_Combination), linewidth=1.2) +
        ggplot2::scale_color_manual(
            values=c(`TRUE`="black", `FALSE`="#E0E0E0")) +
        ggplot2::guides(color="none", fill="none") +
        .upsetTheme + ggplot2::xlab('') + ggplot2::ylab('') +
        ggplot2::theme(
            panel.background=ggplot2::element_blank(),
            axis.text.y=ggplot2::element_text(color='black', size=fontSize))
}

.upsetTheme <- ggplot2::theme(
    panel.background=ggplot2::element_rect(fill='white', color='black'),
    axis.text.x=ggplot2::element_blank(),
    axis.ticks.y=ggplot2::element_blank(),
    axis.text.y=ggplot2::element_text(colour='black'),
    axis.ticks.length=grid::unit(0, "pt"),
    axis.title.y=ggplot2::element_blank(),
    axis.title.x=ggplot2::element_blank(),
    axis.line=ggplot2::element_blank(),
    panel.border=ggplot2::element_blank())

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
