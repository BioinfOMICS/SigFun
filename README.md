# SigFun

## Installation

```r
# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install SigFun with vignettes and dependencies
devtools::install_github(
    "BioinfOMICS/SigFun", 
    build_vignettes = TRUE, 
    dependencies = TRUE)
```