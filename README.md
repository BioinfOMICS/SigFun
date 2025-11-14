# SigFun

[![GitHub issues](https://img.shields.io/github/issues/BioinfOMICS/SigFun)](https://github.com/BioinfOMICS/SigFun/issues) [![GitHub license](https://img.shields.io/github/license/BioinfOMICS/SigFun)](https://github.com/BioinfOMICS/SigFun/blob/main/LICENSE)

## Overview

**SigFun** is an innovative computational framework that leverages whole-transcriptome data to systematically analyze multi-gene signature functions at the system level. SigFun overcomes the limitations of traditional single-gene interpretation and small gene set approaches by providing comprehensive functional insights into clinically validated yet mechanistically opaque signatures. SigFun effectively bridges the critical gap between established clinical utility and biological understanding, enabling researchers to uncover the mechanistic basis underlying signature predictive power in precision medicine applications.

### Key Features

-   **System-level analysis**: Utilizes whole-transcriptome data for comprehensive functional interpretation
-   **Flexible input**: Supports both numeric and binary signature formats
-   **Rich visualization**: 12+ visualization functions for multi-perspective result interpretation
-   **Streamlined workflow**: One-click analysis from SummarizedExperiment input to biological insights

## Vignettes at a Glance

-   **QuickStart with Streamlined Workflow** - Minimal working example *and* a demonstration of the streamlined `sig2Fun` workflow using built-in test data. *Start here to verify installation and understand the basic workflow.*

-   **Data Preparation** - How to build the required `SummarizedExperiment`: expression, `rowData`, `colData`, and `t2g` (with quick validators). *Start here if your data isn’t in SE format yet.*

-   **Stepwise Workflow (`sigCor` → GSEA → plots)** - Runs correlation and enrichment separately so you can inspect/modify `cor.df`, change ranking metrics, or swap enrichment settings. *Use for advanced customization.*

-   **Visualization Functions** - What each SigFun plotting function does (bar, dot, heat, cnet, emap, tree, ridge, lollipop, UpSet, gsea, chord diagram) with concise usage examples. *Use when crafting publication figures.*

-   **Custom Gene Ranking** - Plug in your own gene-level stats (e.g., `log2FC`, `zscore`) by attaching a `cor.df` and setting `ranking.method`. *Use when you already computed rankings externally.*

## Installation

### Prerequisites

Ensure you have R 4.4.0 or later installed. You can download R from <http://www.r-project.org>. If you are a Windows user, install RTools following the guide in the CRAN official wed site <https://cran.r-project.org/bin/windows/>.

### Install Dependencies

#### Install Bioconductor packages

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "org.Hs.eg.db",
    "BiocManager",
    "enrichplot", 
    "clusterProfiler",
    "S4Vectors",
    "DOSE",
    "SummarizedExperiment",
    "BiocStyle"
))
```

#### Install CRAN packages

``` r
install.packages(c(
    "msigdbr",
    "ggridges",
    "DT",
    "pandoc",
    "dplyr",
    "ggplot2", 
    "ggpubr",
    "tibble",
    "tidyr",
    "scales",
    "stringr",
    "GseaVis",
    "circlize",
    "cli",
    "forcats",
    "igraph",
    "randomcoloR",
    "yulab.utils",
    "devtools",
    "roxygen2",
    "rmarkdown",
    "testthat",
    "knitr",
    "tidyverse"
))
```

### Install SigFun

``` r
# Update repositories
options(repos = c(
    CRAN = "https://cloud.r-project.org/",
    BiocManager::repositories()))

# Install SigFun from GitHub
devtools::install_github(
    "BioinfOMICS/SigFun", 
    build_vignettes = TRUE, 
    dependencies = TRUE)
```

## Quick Start

### Load Required Libraries

``` r
library(SigFun)
library(dplyr)
```

### Load Demo Dataset

``` r
# Load example data
data("demo_GSE181574")

# Construct SummarizedExperiment object
GSE181574.sigfun <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = as.matrix(expr.data)),
    rowData = S4Vectors::DataFrame(mapping, row.names = mapping$ensg_id),
    colData = S4Vectors::DataFrame(SIG_MAT))

# Check the structure
show(GSE181574.sigfun)
```

### Run SigFun Analysis

``` r
# Perform signature functional analysis
# Note: This demo uses a binary (1/0) classification signature
GSE181574.sigfun.res <- sig2Fun(
    GSE181574.sigfun,
    cor.method = "logit",  # Use "logit" for binary signatures
    t2g = t2g,             # Ontology dataset
    strings = c("GOBP", "GOCC", "GOMF", "KEGG", 
               "REACTOME", "WP", "HALLMARK", "SIGNALING"))

# View results summary
show(GSE181574.sigfun.res)
```

### View Enrichment Results

``` r
# Extract GSEA results
GSEA_result <- GSE181574.sigfun.res@metadata$gseaResult@result

# View top positively enriched GO Biological Process terms
GSEA_result %>% 
    dplyr::slice(grep("GOBP", GSEA_result$ID)) %>%
    dplyr::arrange(desc(NES)) %>%
    dplyr::select(ID, Description, NES, pvalue, p.adjust, qvalue) %>%
    head(10)

# View top negatively enriched terms  
GSEA_result %>% 
    dplyr::slice(grep("GOBP", GSEA_result$ID)) %>%
    dplyr::arrange(NES) %>%
    dplyr::select(ID, Description, NES, pvalue, p.adjust, qvalue) %>%
    head(10)
```

### Visualization

``` r
# Generate heatmap visualizations
All_heatmaps <- GSE181574.sigfun.res@metadata$heatmap

# Display GO Biological Process heatmap
All_heatmaps$GOBP

# Display HALLMARK pathway heatmap  
All_heatmaps$HALLMARK
```

The visualization includes: - **Carplot (left)**: NES value visualization - **NES (middle-left)**: Statistical results with NES and p-values - **Enrichment plot (middle-right)**: Gene distribution across the transcriptome - **Function name (right)**: Official pathway/function names

### Additional Visualization Functions

SigFun provides 12+ visualization functions for comprehensive result interpretation:

1.  `plot_heat()` — Integrative multi-panel heatmap summarizing pathways by NES, p-values, and enrichment curves.

2.  `barPlot()` — Bar chart ranking top enriched pathways by significance or enrichment score.

3.  `chordPlot()` — Circular diagram linking genes to enriched pathways to show shared and unique memberships.

4.  `cnetPlot()` — Category–gene network visualizing shared genes and connectivity among pathways.

5.  `dotPlot()` — Dot chart encoding significance and gene ratio across pathways.

6.  `emapPlot()` — Enrichment map network showing similarity between pathways based on shared genes.

7.  `gseaPlot()` — Per-pathway running enrichment score curve with ranked gene metrics.

8.  `heatPlot()` — Simple gene × pathway dot heatmap encoding correlation and significance.

9.  `lollipopPlot()` — NES-oriented lollipop chart comparing enrichment magnitude and direction.

10. `ridgePlot()` — Ridge density plot showing gene ranking distributions within each pathway.

11. `treePlot()` — Hierarchical clustering tree grouping pathways by gene overlap or similarity.

12. `upsetPlot()` — UpSet diagram illustrating intersections of genes among enriched pathways.

## Method Selection

Choose the appropriate correlation method based on your signature type:

-   **Numeric signatures**: Use `cor.method = "spearman"` (default), `"pearson"`, or `"kendall"`
-   **Binary signatures**: Use `cor.method = "logit"` (univariate logistic regression)

## Input Data Requirements

SigFun requires a `SummarizedExperiment` object with:

1.  **Assay**: Gene expression matrix (genes × samples)
2.  **RowData**: Gene information including `ensg_id`, `gene_symbol`, and `gene_biotype`
3.  **ColData**: Sample information with `sample_id` and signature `value`

## License

This project is licensed under the [MIT License](LICENSE).
