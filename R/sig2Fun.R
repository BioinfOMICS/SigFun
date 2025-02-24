#' @title sig2Fun
#' @description This function is explore the functions related to signature by a
#' surrogate of whole transcriptome data.
#' @param SE_data Summarized Experiment Object. Input contains "colData",
#' "assays", and "elementMetadata".
#' Use SummarizedExperiment::SummarizedExperiment to integrate your data.
#' The example code of constructing the input object:
#' SE_data <- SummarizedExperiment::SummarizedExperiment(
#' assays=list(abundance=as.matrix(assay.data)),
#' rowData=S4Vectors::DataFrame(rowData.data, row.names=rowData.data$sample_id),
#' colData=colData.data)
#' assay.data is the expression value where gene in column and sample in row
#' rowData is a vector of sample ID.
#' colData.data is the gene information containing three columns:
#' "ensg_id" (ensemble ID), "gene_symbol" (official gene symbol), and
#' "gene_biotype" (annotate the gene is "protein_coding" or not).
#' Use data(demo) to see the example.
#' @param ranking.method Character. "stat" (default) or "pval"
#' @param species Character. "human" (default) or "mouse"
#' @param cor.method Character. "spearman" (default), "pearson", or "kendall"
#' @param pathways_all list file. Loaded by fgsea::gmtPathways. Should be
#' obtained from MSigDB, such as msigdb.v2023.1.Hs.symbols.gmt
#' @param output_path Character. Location of the output path.
#' @param topN Integer. Number of significant functions in hte output bar plots.
#' Default is 10.
#' @param Z.transform logic. Do z-transform on GENE_MAT or not (default: FALSE).
#' @param significat_type Character. Filter or statistical significance.
#' "pval" for p-value and "qval" for Q-value
#' @param strings vector. Pathways header want to be plotted.
#' Default c("GOBP", "GOCC", "GOMF", "KEGG", "REACTOME", "WP") for
#' msigdb.v2023.1.Hs.symbols.gmt
#' @param plot_out logic. Plot the heatmap to output_path (Default: FALSE).
#' @return Return a SummarizedExperiment object containing analysis results.
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom stats cor
#' @importFrom fgsea gmtPathways
#' @importFrom methods new
#' @importFrom grDevices pdf dev.off
#' @importFrom utils write.csv
#' @include utility.R
#' @include sigCor.R
#' @include sig2GSEA.R
#' @include plot_bar.R
#' @include plot_heat.R
#' @export
#' @examples
#' data("demo")
#' output_path <- tempdir()
#' res <- sig2Fun(SE_data, ranking.method="stat", species="human",
#' cor.method="spearman", pathways.all, output_path,
#' topN=10, Z.transform=FALSE, significat_type="pval",
#' strings=c("KEGG"))


sig2Fun <- function(SE_data, ranking.method="stat", species="human",
                    cor.method="spearman", pathways_all, output_path,
                    topN=10, Z.transform=FALSE, significat_type="pval",
                    strings=c("GOBP","GOCC","GOMF","KEGG","REACTOME","WP"),
                    plot_out=FALSE) {
    .classCheck(SE_data, "SummarizedExperiment")
    .typeCheck(output_path, "character")
    .classCheck(pathways_all, "list")

    if (!dir.exists(output_path)) {
        dir.create(output_path, recursive=TRUE, showWarnings =FALSE)
    }
    Sys.chmod("my_folder", mode="755", use_umask=TRUE)
    SE_data.cor <- SE_data
        if(!("cor.df" %in% names(SE_data@metadata))){
            SE_data.cor <- sigCor(SE_data=SE_data, cor.method=cor.method,
            output_path=output_path, Z.transform=Z.transform)
            metadata(SE_data) <- list(cor.df=metadata(SE_data.cor)$cor.df)
        }

        if(!("fgseaRes" %in% names(SE_data@metadata))){
            SE_data.fgsea <- sig2GSEA(SE_data.cor=SE_data.cor,
                ranking.method=ranking.method, output_path=output_path,
                pathways.all=pathways_all)
            .summary_gsea(metadata(SE_data.fgsea)$fgsea,
            ranking.method=ranking.method, output_path=output_path)
            metadata(SE_data) <- list(fgsea=metadata(SE_data.fgsea)$fgsea)
        }

    if(plot_out){
    barplots <- plot_bar(SE_data.fgsea=SE_data, output_path=output_path,
    topN=topN, significat_type=significat_type, strings=strings)
    metadata(SE_data) <- list(barplots=barplots)

    heatmap <- plot_heat(SE_data.fgsea=SE_data, output_path=output_path,
        strings=strings, significat_type=significat_type, topN=topN,
        pathways.all=pathways_all, ranking.method=ranking.method)
    metadata(SE_data) <- list(heaymap=heatmap)
    }

    return(SE_data)
}
