#' @title SigFun
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
#' @param pathways.all list file. Loaded by fgsea::gmtPathways. Should be
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
#' @export
#' @examples
#' data("demo")
#' output_path <- tempdir()
#' res <- SigFun(SE_data, ranking.method = "stat", species="human",
#' cor.method = "spearman", pathways.all, output_path,
#' topN=10, Z.transform=FALSE, significat_type="pval",
#' strings=c("KEGG"))


SigFun <- function(SE_data, ranking.method = "stat", species="human",
                    cor.method = "spearman", pathways.all, output_path,
                    topN=10, Z.transform=FALSE, significat_type="pval",
                    strings=c("GOBP","GOCC","GOMF","KEGG","REACTOME","WP"),
                    plot_out=FALSE) {

    # Step 0. setup environment: create output folder ####
    if (!dir.exists(output_path)) {
        dir.create(output_path, recursive = TRUE, showWarnings =FALSE)
    }
    system2(paste0("chmod 775 -R ", output_path))
    res.all <- SE_data
    # Step 1. SigCor: run gene expression correlation with signature patterns ##

        if(!("cor.df" %in% names(SE_data@metadata))){
            SE_data.cor <- SigCor(SE_data=SE_data, cor.method=cor.method,
            output_path=output_path, Z.transform=Z.transform)
            res.all@metadata$cor.df <- SE_data.cor@metadata$cor.df
        }

    # Step 2. Sig2GSEA: Run fGSEA using gene ranking values load corRES ####
        if(!("fgseaRes" %in% names(SE_data@metadata))){
            SE_data.fgsea <- Sig2GSEA(SE_data.cor=SE_data.cor,
                ranking.method=ranking.method, output_path=output_path,
                pathways.all=pathways.all)
        # write result tables: GSEA_pathway_stat_pvalue.txt; qvalue.txt; all.txt
            .summary_gsea(SE_data.fgsea@metadata$fgseaRes,
            ranking.method=ranking.method, output_path=output_path)
            res.all@metadata$fgseaRes <- SE_data.fgsea@metadata$fgseaRes
        }

    # Step 3. plot.bar ####
    if(plot_out){
    barplots <- plot_bar(SE_data.fgsea=res.all, output_path=output_path,
    topN=topN, significat_type=significat_type, strings = strings)
    res.all[["barplots"]] <- barplots

    # Step 4. plot summary heatmap table  ####
    heatmap <- plot_heat(SE_data.fgsea=res.all, output_path=output_path,
        strings=strings, significat_type=significat_type, topN=topN,
        pathways.all=pathways.all, ranking.method=ranking.method)
    res.all[["heatmap"]] <- heatmap
    }

    return(res.all)
}
