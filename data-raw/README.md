# Regenerating packaged vignette data

The script in this directory rebuilds the CTSS example objects used by the vignettes:

- `inst/extdata/ctss.rds`
- `inst/extdata/ctss_clean.rds`

Run it from the package root with:

```sh
make extdata
```

or directly with:

```sh
Rscript data-raw/make-ctss-extdata.R
```

The script uses the bundled `inst/extdata/design_matrix_first10pct.tsv` file and BigWig files in `inst/extdata/cage_bw/`, following the same processing steps shown in `vignettes/ctss-processing.Rmd`.
