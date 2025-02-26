# SigFun
SigFun is an R package designed to streamline the analysis of transcriptomic data in relation to specific gene signatures. 
It provides an automated workflow for analyzing gene signatures and their functional implications in transcriptomic datasets. 
The package integrates Gene Set Enrichment Analysis (GSEA) with visualization tools to help researchers understand the biological 
pathways and processes associated with their gene signatures of interest.

## Installation 
Before the installation, you should first install the following packages from Bioconductor:

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("SigFun")
``` 

After conducting the above step, now you can load in our package and start 
using it!
