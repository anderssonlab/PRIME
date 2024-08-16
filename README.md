# PRIME '

## `PRIME` - regulatory element analysis using transcription start site data

`PRIME` gathers a suite of `R` functions for:
* analyzing transcription start site (TSS) data, e.g. CAGE, or other RNA 5' end data like GRO/PRO-cap
* analyzing divergent transcription
* calling enhancers from TSS data
* subsampling of data, e.g. for saturation analyses
* decomposition of TSS clusters for core promoter analysis
* calculation of genomic background (noise) expression
* normalization across libraries versus GC content
* calculation of bias in expression coverage and expression support versus batch

## Installation:

Install `PRIME` directly from GitHub using `devtools`:

```
devtools::install_github("anderssonlab/PRIME")
```

You may need to manually install dependencies (See `DESRIPTION`).

## Notes:
* PRIME uses [`CAGEfightR`](https://github.com/MalteThodberg/CAGEfightR) as a base for data loading and handling
* R package `bcp` is not vailable from CRAN for `R` >= 4.3. To install `bcp`, use `install_github("swang87/bcp")`

## Contributors
* Robin Andersson
* Natsuda Navamajiti
* Hjorleifur Einarsson
* Robert Krautz
* Nicolas Alcaraz

## TO DO
* write vignettes

## Future development
* ATAC-CAGE integration
* QC tools and visualization
