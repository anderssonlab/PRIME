# PRIME

[![R-CMD-check](https://github.com/anderssonlab/PRIME/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/anderssonlab/PRIME/actions/workflows/R-CMD-check.yaml)

## Quantitative analysis of transcription initiation data

**PRIME** is an R package for **quantitative analysis of transcription initiation data** from **CAGE** and related **TSS/CTSS assays**. It extends **[CAGEfightR](https://github.com/MalteThodberg/CAGEfightR)** with tools for **regulatory element characterization**, including transcription initiation **complexity**, **signal saturation**, **background noise**, **divergent transcription**, and **core promoter architecture**. PRIME also provides an interface to a **LightGBM model** for data-driven regulatory element scoring (via **[PRIMEmodel](https://github.com/anderssonlab/PRIMEmodel)**).

PRIME is designed to work alongside **CAGEfightR** and uses the same Bioconductor data structures (**SummarizedExperiment** / **RangedSummarizedExperiment** and **GRanges**). PRIME adds additional quantitative and modeling utilities while staying fully compatible with Bioconductor genomic infrastructure.

## Core functionality

### CTSS-level quantification & processing
- CTSS-level quantification (via `CAGEfightR::quantifyCTSSs()`)
- expression summary, normalization, subsampling
- transcription initiation metrics
  - complexity (e.g., dispersion/entropy-style summaries depending on analysis)
  - saturation / downsampling-based analyses
  - background noise estimation

### Regulatory element characterization
- divergent transcription detection and quantification
- tag cluster analysis utilities (decomposition and downstream quantification)
- strand balance and other profile-derived summaries (depending on workflow)

### Core promoter analysis
- promoter decomposition
- positional dispersion summaries (e.g., width/dispersion-style measures)
- initiator sequence patterns (INR-like classifications)

### Profile-based analysis
- signal aggregation around genomic features
- window-based summarization and heatmap-style matrices

### Machine learning interface
- interfaces to score candidate regulatory elements using the **PRIME LightGBM model** (distributed via **[PRIMEmodel](https://github.com/anderssonlab/PRIMEmodel)**)

## PRIME toolkit

The **PRIME toolkit** consists of three interconnected tools for the analysis of transcription initiation data (e.g., CAGE):

| Tool | Type | Purpose |
|---|---|---|
| [PRIMEprep](https://github.com/anderssonlab/PRIMEprep) | Bash pipeline | Raw FASTQ → QC → trimming → mapping → BigWig |
| [PRIME](https://github.com/anderssonlab/PRIME) | R package | CTSS quantification, divergent loci, promoter decomposition, normalization, noise estimation |
| [PRIMEmodel](https://github.com/anderssonlab/PRIMEmodel) | R package + Python | Genome-wide prediction / scoring of regulatory elements |

## Input data

PRIME works with standard Bioconductor genomic data structures:

- **CTSS-level data**: typically a `RangedSummarizedExperiment` produced by CAGEfightR (row ranges are CTSS positions; assays contain counts/TPM).
- **Tag clusters / loci / regions**: typically `GRanges` (or `RangedSummarizedExperiment` objects where `rowRanges()` are regions).
- PRIME functions are generally compatible with `SummarizedExperiment` + `GRanges` workflows and integrate with the Bioconductor ecosystem.

## Getting started

PRIME is a toolbox (not a single rigid pipeline). The fastest way to get started is to follow the vignettes on the PRIME website:

- Articles index: https://anderssonlab.org/PRIME/articles/
- Installing PRIME: https://anderssonlab.org/PRIME/articles/installation.html
- CTSS processing & QC: https://anderssonlab.org/PRIME/articles/ctss-processing.html
- Tag cluster decomposition: https://anderssonlab.org/PRIME/articles/tag-cluster-decomposition.html
- Divergent loci: https://anderssonlab.org/PRIME/articles/divergent-loci.html
- Normalization & batches: https://anderssonlab.org/PRIME/articles/normalization-batches.html
- Noise estimation: https://anderssonlab.org/PRIME/articles/noise-estimation.html
- Regulatory element prediction (PRIMEmodel): https://anderssonlab.org/PRIME/articles/prediction.html
- End-to-end workflow: https://anderssonlab.org/PRIME/articles/end-to-end-workflow.html

## Relationship to CAGEfightR

PRIME is designed to work alongside **CAGEfightR**. It uses the same data structures (`SummarizedExperiment` and `GRanges`) and extends CAGEfightR with additional quantitative and modeling utilities for regulatory element analysis.

CAGEfightR repository: https://github.com/MalteThodberg/CAGEfightR

## PRIME model (PRIMEmodel)

- **PRIME (this package)** is the analysis toolkit (CTSS QC, quantification helpers, complexity/noise utilities, tag cluster and promoter analysis helpers).
- **PRIMEmodel** distributes the trained LightGBM model and provides genome-wide (or focal) scoring of candidate regulatory elements from CTSS profiles.

PRIMEmodel repository: https://github.com/anderssonlab/PRIMEmodel  
PRIMEmodel website: https://anderssonlab.org/PRIMEmodel/

## Example applications

- Compare transcription initiation complexity across conditions
- Identify and characterize divergently transcribed loci
- Analyze core promoter architecture via decomposition and dispersion measures
- Perform profile-based analyses around genomic features (heatmap-style and window summaries)
- Score candidate regulatory elements using the PRIME model (via PRIMEmodel)

## Installation

PRIME depends on a mix of CRAN, Bioconductor, and GitHub packages. For detailed instructions, see:

- https://anderssonlab.org/PRIME/articles/installation.html

## Documentation

- Website: https://anderssonlab.org/PRIME/
- Articles (vignettes): https://anderssonlab.org/PRIME/articles/
- Reference: https://anderssonlab.org/PRIME/reference/

## Developers

- [Robin Andersson](https://github.com/anderssonrobin)
- [Natsuda Navamajiti](https://github.com/natsnav)
- [Hjorleifur Einarsson](https://github.com/HjolliEin)
- [Robert Krautz](https://github.com/robertkrautz)
- [Nicolas Alcaraz](https://github.com/satroz)
