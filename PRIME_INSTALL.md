# PRIME Installation Preparation Guide

This guide explains how to set up the required **R** environments to use the `PRIME` R package.
The R package `PRIME` depends on a mix of CRAN and Bioconductor, and GitHub packages.

```bash
R
```

1. Install required CRAN packages
```r
install.packages(c(
  "assertthat",
  "data.table",
  "argparse",
  "arrow",
  "Matrix",
  "caTools",
  "igraph",
  "zoo"
))
```

2. Install BiocManager (if not already installed)
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
```

3. Install required Bioconductor packages
```r
## General principle: Installing `CAGEfightR` with:

BiocManager::install("CAGEfightR")

## will automatically pull in dependencies.
## You do NOT need to install these manually unless an error occurs.**
```
If errors occur, install in layers:
```r
## core:
BiocManager::install("S4Vectors", "IRanges", "GenomeInfoDb")

## data structures:
BiocManager::install("SummarizedExperiment", "GenomicRanges")

## utilities
BiocManager::install("BiocParallel", "BSgenome", "rtracklayer", "Hmisc")

## CAGEfightR
BiocManager::install("CAGEfightR")
```

4. Install devtools
```r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
```

5. Install bcp from github
```r
devtools::install_github("swang87/bcp")
```

6. Install PRIME
```r
devtools::install_github("anderssonlab/PRIME")
```
---

© 2025 PRIME setup protocol
