# PRIME '

## `PRIME` - regulatory element analysis using transcription start site data

`PRIME` gathers a suite of `R` functions for:
* analyzing transcription start site (TSS) data, e.g. CAGE, or other RNA 5' end data like GRO/PRO-cap
* analyzing divergent transcription
* calling enhancers and promoters from TSS data
* subsampling of data, e.g. for saturation analyses
* decomposition of TSS clusters for core promoter analysis
* calculation of genomic background (noise) expression
* normalization across libraries
* dealing with experimental batches

## Installation

The `R` package `PRIME` depends on a mix of CRAN, Bioconductor, and GitHub packages. Installation also requires a specific `python` environment setup. Please refer to [Installation instructions](https://github.com/anderssonlab/PRIME/blob/main/INSTALL_PRIME.md) for details.

## Associated repos

* [PRIMEprep](https://github.com/anderssonlab/PRIMEprep) Preprocessing and mapping of CAGE sequencing data. Recommended to prepare data for `PRIME` analysis.
* [PRIMEloci](https://github.com/anderssonlab/PRIMEloci) Prediction of regulatory elements from transcription initiation data. PRIMEloci is integrated in the `PRIME` R package.

## Contributors
* [Robin Andersson](https://github.com/anderssonrobin)
* [Natsuda Navamajiti](https://github.com/natsnav)
* [Hjorleifur Einarsson](https://github.com/HjolliEin)
* [Robert Krautz](https://github.com/robertkrautz)
* [Nicolas Alcaraz](https://github.com/satroz)
