#' @title sigCor
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
#' @param cor.method Character. "spearman" (default), "pearson", or "kendall"
#' @param Z.transform logic. Do z-transform on GENE_MAT or not (default: FALSE).
#' @param output_path Character. Location of the output path.
#' @return Return a SummarizedExperiment object containing analysis results.
#' SigCor will add "cor.df" to the metadata of the input object which is the
#' necessary input for Sig2GSEA.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assay
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather pivot_wider
#' @importFrom dplyr filter select distinct
#' @importFrom stats cor
#' @importFrom readr write_delim
#' @importFrom utils write.csv
#' @include utility.R
#' @export
#' @examples
#' data("demo")
#' output_path <- tempdir()
#' SE_data.cor <- sigCor(
#' SE_data, cor.method="spearman", output_path=output_path,
#' Z.transform=FALSE)


sigCor <- function(SE_data, output_path,
                    cor.method="spearman", Z.transform=FALSE) {
    signature.obj <- as.data.frame(SummarizedExperiment::rowData(SE_data))

    RNA_exp_sql <- as.data.frame(SummarizedExperiment::assay(SE_data)) %>%
        tibble::rownames_to_column("sample_id") %>%
        tidyr::gather(-sample_id, key="ensg_id", value="value")
    corRES.path <- file.path(output_path, "corRES.txt")
    S_ID <- intersect(signature.obj$sample_id, unique(RNA_exp_sql$sample_id))
    RNA_exp_sql.sid <- RNA_exp_sql %>% dplyr::filter(sample_id %in% S_ID)
    signature.obj <- signature.obj %>% dplyr::filter(sample_id %in% S_ID)
    RNA_exp_sql.FPKM_UQ.tmp <- RNA_exp_sql.sid %>%
        dplyr::select(sample_id, ensg=ensg_id, value)
    tmp.data.wide <- RNA_exp_sql.FPKM_UQ.tmp %>%
        dplyr::filter(ensg != "") %>%
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
    } else {
        apply(pattern.genes, 2, .z_score_cal)
    }
    } else {
        pattern.genes.norm <- if (length(rm.index) > 0) {
        pattern.genes[, -rm.index]
        } else {
            apply(pattern.genes, 2, .z_score_cal)
        }
    }
    cor.list <- apply(pattern.genes.norm, 2, function(x) {
        cor(x, pattern.signature,
        method=cor.method, use="pairwise.complete.obs")
    })
    cor.df <- data.frame(gene=names(cor.list), cor=as.numeric(cor.list))
    readr::write_delim(cor.df, corRES.path, delim="\t")
    SE_data@metadata$cor.df <- cor.df
    return(SE_data)
}
