````markdown name=INSTALL.md
# PRIME Installation Guide

This guide describes how to install the **PRIME** R package and its dependencies.

PRIME depends on a mix of **CRAN**, **Bioconductor**, and **GitHub** packages.

---

## Quick start (recommended)

Start an R session:

```bash
R
```

### 1) Install BiocManager

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
```

### 2) Install Bioconductor dependencies (recommended: start with CAGEfightR)

PRIME imports `CAGEfightR`. Installing it will pull in most Bioconductor dependencies automatically.

```r
BiocManager::install("CAGEfightR", ask = FALSE, update = FALSE)
```

### 3) If installation fails due to `pryr` (archived on CRAN)

Some dependency stacks require `pryr (>= 0.1.3)`. `pryr` was archived on CRAN on **2026-01-30**.

If you see an error mentioning `pryr`, install it from the CRAN archive:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_url(
  "https://cran.r-project.org/src/contrib/Archive/pryr/pryr_0.1.6.tar.gz",
  upgrade = "never"
)
```

Then retry installing `CAGEfightR`:

```r
BiocManager::install("CAGEfightR", ask = FALSE, update = FALSE)
```

### 4) Install devtools (or remotes)

You can install PRIME using either `devtools` or `remotes`. If you already installed `remotes` above, you can skip installing `devtools`.

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
```

### 5) Install `bcp` from GitHub

`bcp` is installed from GitHub:

```r
devtools::install_github("swang87/bcp")
```

### 6) Install PRIME from GitHub

```r
devtools::install_github("anderssonlab/PRIME")
```

---

## Optional: install a small set of commonly-used CRAN packages

Most CRAN dependencies for PRIME will be installed automatically when you install PRIME. If you want to proactively install a few common packages (helpful on some systems), you can do:

```r
install.packages(c(
  "assertthat",
  "Matrix",
  "data.table",
  "caTools",
  "dplyr",
  "foreach",
  "ggplot2",
  "igraph",
  "magrittr",
  "purrr",
  "stringr",
  "tibble",
  "tidyr",
  "tidyselect",
  "zoo"
))
```

---

## Notes on R / Bioconductor versions

- PRIME should be installable on **newer versions of R**. In general, you should use the **Bioconductor release that matches your R version** (BiocManager will normally pick the correct one automatically).
- The PRIME GitHub Actions workflow used to build the pkgdown site is pinned to:
  - **R 4.2.2**
  - **Bioconductor 3.16**

If you need to reproduce that exact environment (e.g., for server installs or CI parity), you can explicitly pin Bioconductor like this:

```r
BiocManager::install(version = "3.16", ask = FALSE, update = FALSE)
```

If you are using a newer R version, do **not** force Bioconductor 3.16; instead, install the Bioconductor version recommended for your R version.

---

2026 PRIME installation guide
````