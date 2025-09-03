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
