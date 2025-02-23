% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig2GSEA.R
\name{sig2GSEA}
\alias{sig2GSEA}
\title{sig2GSEA}
\usage{
sig2GSEA(SE_data.cor, ranking.method, output_path, pathways.all, NPROC = 1)
}
\arguments{
\item{SE_data.cor}{Summarized Experiment Object. In the metadata of the input
, it must be a Summary Experiment Object and should contain the following
objects "cor.df". "cor.df" is the result of SigCor.
Add the cor.df to the metadata by running SigCor.}

\item{ranking.method}{Character. "stat" (default) or "pval"}

\item{output_path}{Character. Location of the output path.}

\item{pathways.all}{gmt file. loading by fgsea::gmtPathways}

\item{NPROC}{Parallelization parameter used in bplapply (default=1).
Integer Number of permutations to do (default=1000).}
}
\value{
Return a SummarizedExperiment object containing analysis results.
SigCor will add "fgseaRes" to the metadata of the input object which is the
necessary input for plot_bar and plot_heat.
}
\description{
This is a key component, suggesting the package performs Gene
Set Enrichment Analysis (GSEA). It uses the SigCor results (correlation
statistics) and/or p-values, along with the ontology profile, to identify
enriched pathways or gene sets associated with the input signature.
}
\examples{
data(sig2GSEA_input)
temp_dir <- tempdir()
result <- sig2GSEA(SE_data.cor=test_SE, ranking.method="stat",
output_path=temp_dir, pathways.all=test_pathways, NPROC=2)
}
