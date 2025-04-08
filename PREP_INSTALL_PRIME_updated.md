# 🧬 PRIME Installation Preparation Guide

This guide explains how to set up the required **R** and **Python** environments to use the `PRIME` R package. The Python environment is also compatible with `PRIMEloci()` if used.

---

## 📦 R Installation for PRIME

The R package `PRIME` depends on a mix of CRAN and Bioconductor packages, and a custom GitHub version of `bcp`.

### ✅ Full R Setup

```r
# 1. (macOS only) Optional: Fix Fortran flags for compiling 'bcp' from GitHub
makevars_path <- path.expand("~/.R/Makevars")
if (!file.exists(makevars_path)) {
  dir.create(dirname(makevars_path), showWarnings = FALSE)
  cat(
    "FC = /usr/local/bin/gfortran\n",
    "F77 = /usr/local/bin/gfortran\n",
    "FLIBS = -L/usr/local/lib/gcc/13 -lgfortran -lquadmath\n",
    file = makevars_path
  )
} else {
  message("~/.R/Makevars already exists. Please verify Fortran flags.")
}
```

```r
# 2. Install required CRAN packages
install.packages(c(
  "future",
  "future.apply",
  "foreach",
  "argparse",
  "doParallel",
  "reticulate",
  "arrow",
  "Matrix",
  "tools"
))

# 3. Install BiocManager (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 4. Install required Bioconductor packages
BiocManager::install(c(
  "CAGEfightR",
  "SummarizedExperiment",
  "S4Vectors",
  "GenomicRanges",
  "BiocParallel",
  "IRanges",
  "GenomeInfoDb",
  "BSgenome",
  "rtracklayer",
  "sparseMatrixStats",
  "DAPAR",
  "DAPARdata"
))

# 5. Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# 6. Install bcp from GitHub (custom version)
devtools::install_github("swang87/bcp")
```

---

## 🐍 Python Environment Setup (Compatible with PRIMEloci)

You can prepare the Python environment in **two ways**:

---

### 🧬 Option 1: Conda (recommended for reproducibility and flexibility)

> 💡 **Note:** You only need to set up one of these environments (NumPy 1.x or NumPy 2.x).  
> Both are fully compatible with the PRIME package and its Python integration.

`environment_np1.yml` and `environment_np2.yml` are included in the `inst/` folder.

#### How to create the environment:

```bash
cd PRIME

# For NumPy 1.x
conda env create -f inst/conda_env/environment_np1.yml

# OR for NumPy 2.x
conda env create -f inst/conda_env/environment_np2.yml

# Activate the environment
conda activate prime-env-numpy1  # or prime-env-numpy2
which python
# Copy this path to use as python_path in R when calling run_PRIMEloci_example() or run_PRIMEloci_facet_example()
```

#### Example in R:

```r
run_PRIMEloci_example(python_path = "~/.conda/envs/prime-env-numpy1/bin/python")
run_PRIMEloci_facet_example(python_path = "~/.conda/envs/prime-env-numpy1/bin/python")
```

---

### 🧪 Option 2: Virtualenv

#### `requirements.txt`

```txt
numpy==1.26.4
scipy==1.11.4
pandas==2.1.4
joblib==1.3.2
scikit-learn==1.3.2
lightgbm==4.1.0
pyarrow==14.0.2
fastparquet==2023.8.0
```

#### Setup instructions

```bash
python3 -m venv ~/prime_env_np1
source ~/prime_env_np1/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
which python
# Copy this path to use as python_path in R when calling run_PRIMEloci_example() or run_PRIMEloci_facet_example()
```

#### Example in R:

```r
run_PRIMEloci_example(python_path = "~/prime_env_np1/bin/python")
run_PRIMEloci_facet_example(python_path = "~/prime_env_np1/bin/python")
```

---

## 🧠 Tips

- Always check the current Python path with `which python` **after** activating your environment.
- Use that full path in the `python_path` argument in R.
- Both conda environments are compatible with `PRIMEloci()` if needed later.

---

© 2025 PRIME setup protocol | Supports PRIMEloci-compatible Python integration
