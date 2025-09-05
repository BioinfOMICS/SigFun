# SigFun

[![R-CMD-check](https://github.com/BioinfOMICS/SigFun/workflows/R-CMD-check/badge.svg)](https://github.com/BioinfOMICS/SigFun/actions)
[![GitHub issues](https://img.shields.io/github/issues/BioinfOMICS/SigFun)](https://github.com/BioinfOMICS/SigFun/issues)
[![GitHub license](https://img.shields.io/github/license/BioinfOMICS/SigFun)](https://github.com/BioinfOMICS/SigFun/blob/main/LICENSE)

## Overview

**SigFun** is an innovative computational framework that leverages whole-transcriptome data to systematically analyze multi-gene signature functions at the system level. SigFun overcomes the limitations of traditional single-gene interpretation and small gene set approaches by providing comprehensive functional insights into clinically validated yet mechanistically opaque signatures. SigFun effectively bridges the critical gap between established clinical utility and biological understanding, enabling researchers to uncover the mechanistic basis underlying signature predictive power in precision medicine applications.

### Key Features

- **System-level analysis**: Utilizes whole-transcriptome data for comprehensive functional interpretation
- **Flexible input**: Supports both numeric and binary signature formats
- **Rich visualization**: 12+ visualization functions for multi-perspective result interpretation
- **Streamlined workflow**: One-click analysis from SummarizedExperiment input to biological insights

## Installation

### Prerequisites

Ensure you have R 4.4.0 or later installed. You can download R from [http://www.r-project.org](http://www.r-project.org).
If you are a Windows user, install RTools following the guide in the CRAN official wed site [https://cran.r-project.org/bin/windows/](https://cran.r-project.org/bin/windows/).

### Install Dependencies

#### Install Bioconductor packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
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

```r
install.packages(c(
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

```r
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

```r
library(SigFun)
library(dplyr)
```

### Load Demo Dataset

```r
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

```r
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

```r
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

```r
# Generate heatmap visualizations
All_heatmaps <- GSE181574.sigfun.res@metadata$heatmap

# Display GO Biological Process heatmap
All_heatmaps$GOBP

# Display HALLMARK pathway heatmap  
All_heatmaps$HALLMARK
```

The visualization includes:
- **Carplot (left)**: NES value visualization
- **NES (middle-left)**: Statistical results with NES and p-values
- **Enrichment plot (middle-right)**: Gene distribution across the transcriptome
- **Function name (right)**: Official pathway/function names

### Additional Visualization Functions

SigFun provides 12+ visualization functions for comprehensive result interpretation:

- `plot_heat()` — Integrative heatmap visualization
- `barPlot()` — Bar plot visualization  
- `chordPlot()` — Chord plot visualization
- `cnetPlot()` — Concept network visualization
- `dotPlot()` — Dot plot visualization
- `emapPlot()` — Enrichment map visualization
- `gseaPlot()` — Enrichment score plot
- `heatPlot()` — Heat plot visualization
- `lollipopPlot()` — Lollipop plot visualization
- `ridgePlot()` — Ridge plot visualization
- `treePlot()` — Tree plot visualization
- `upsetPlot()` — UpSet plot visualization

## Method Selection

Choose the appropriate correlation method based on your signature type:

- **Numeric signatures**: Use `cor.method = "spearman"` (default), `"pearson"`, or `"kendall"`
- **Binary signatures**: Use `cor.method = "logit"` (univariate logistic regression)

## Input Data Requirements

SigFun requires a `SummarizedExperiment` object with:

1. **Assay**: Gene expression matrix (genes × samples)
2. **RowData**: Gene information including `ensg_id`, `gene_symbol`, and `gene_biotype`
3. **ColData**: Sample information with `sample_id` and signature `value`

## License

This project is licensed under the [MIT License](LICENSE).
