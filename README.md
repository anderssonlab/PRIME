# PRIME '

## PRIME - regulatory element analysis using transcription initiation data

PRIME gathers a suite of R functions for:
* calling enhancers from CAGE data
* analyzing divergent transcription
* subsampling of CAGE data, e.g. for saturation analyses
* decomposition of CAGE tag clusters for core promoter analysis
* calculation of genomic background (noise) expression
* normalization across libraries versus GC content
* calculation of bias in expression coverage and expression support versus batch

## Installation:
```
devtools::install_github("anderssonlab/PRIME")
```

## Notes:
* PRIME uses [CAGEfightR](https://github.com/MalteThodberg/CAGEfightR) as a base for CAGE data loading and handling
* R package bcp is not vailable from CRAN for R >= 4.3. To install bcp, use `install_github("swang87/bcp")`

## Contributors
* Robin Andersson
* Natsuda Navamajiti
* Hjorleifur Einarsson
* Robert Krautz
* Nicolas Alcaraz

## TO DO
* write vignettes
