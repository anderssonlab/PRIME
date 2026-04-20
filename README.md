# PRIME

## Regulatory element analysis using transcription start site data

`PRIME` gathers a suite of `R` functions for:

* analyzing transcription start site (TSS) data, e.g. CAGE, or other RNA 5' end data like GRO/PRO-cap
* analyzing divergent transcription
* calling enhancers and promoters from TSS data
* subsampling of data, e.g. for saturation analyses
* decomposition of TSS clusters for core promoter analysis
* calculation of genomic background (noise) expression
* normalization across libraries
* dealing with experimental batches

## PRIME toolkit

The **PRIME toolkit** consists of three interconnected tools for the analysis of
transcription initiation data (e.g., CAGE):

| Tool | Type | Purpose |
|---|---|---|
| [PRIMEprep](https://github.com/anderssonlab/PRIMEprep) | Bash pipeline | Raw FASTQ → QC → trimming → mapping → BigWig |
| [PRIME](https://github.com/anderssonlab/PRIME) | R package | CTSS quantification, divergent loci, promoter decomposition, normalization, noise estimation |
| [PRIMEmodel](https://github.com/anderssonlab/PRIMEmodel) | R package + Python | Genome-wide prediction of regulatory elements |

## Installation

The `R` package `PRIME` depends on a mix of CRAN, Bioconductor, and GitHub packages. Please refer to [Installation instructions](https://anderssonlab.github.io/PRIME/articles/installation.html) for details.

## Developers

* [Robin Andersson](https://github.com/anderssonrobin)
* [Natsuda Navamajiti](https://github.com/natsnav)
* [Hjorleifur Einarsson](https://github.com/HjolliEin)
* [Robert Krautz](https://github.com/robertkrautz)
* [Nicolas Alcaraz](https://github.com/satroz)
