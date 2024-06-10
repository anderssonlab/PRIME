# PRIME - regulatory element analysis using transcription initiation data

Using CAGEfightR as a base, PRIME gathers a suite of functions for

* calling enhancers from CAGE data
* analyzing divergent transcription
* subsampling of CAGE data, e.g. for saturation analyses
* decomposition of CAGE tag clusters for core promoter analysis
* calculation of genomic background (noise) expression
* normalization across libraries versus GC content
* calculation of bias in expression coverage and expression support versus batch

### Installation:
```
devtools::install_github("anderssonlab/PRIME")
```

Note: to install bcp, use install_github("swang87/bcp")
