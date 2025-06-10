# 🧪 PRIME Installation Preparation Guide

This guide explains how to set up the required **R** and **Python** environments to use the `PRIME` R package. The Python environment is also compatible with `PRIMEloci()` if used.

---

## 📦 R Installation for PRIME

The R package `PRIME` depends on a mix of CRAN and Bioconductor packages, and a custom GitHub version of `bcp`.

### ✅ Full R Setup

```r
# 1. Install required CRAN packages
install.packages(c(
  "future",
  "future.apply",
  "foreach",
  "argparse",
  "doParallel",
  "reticulate",
  "arrow",
  "Matrix",
  "caTools"
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
  "sparseMatrixStats"
))

# 5. Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# 6. Install bcp from GitHub (custom version)
devtools::install_github("swang87/bcp")

# *** specify << lib = "/path/to/R_packakges/" >> in install.packages() and BiocManager::install() if necessary
# *** using remotes package instead of devtools if not successfully install devtools
```

---

## 🐍 Python Environment Setup (Compatible with PRIMEloci functions)

You can prepare the Python environment in **four ways**. The recommended default is the `reticulate`-managed virtualenv using **NumPy 2.x** (use NumPy 1.x only if required for legacy compatibility)

---

### 🌐 Option 1 (Default): Virtualenv via `reticulate` in R

This is the **default and recommended** method for setting up the Python environment using only R. It ensures full compatibility with `PRIME` and `PRIMEloci()` 

#### Setup instructions in R:

```bash
cd PRIME
R
```

```r
library(reticulate)

# Define the environment name
env_name <- "prime-env"

# Create the virtual environment if it doesn't exist
virtualenv_create(envname = env_name)

# Use and configure it
use_virtualenv(env_name, required = TRUE)

# Install Python packages
required_pkgs <- readLines("/PATH/to/PRIME/inst/envfile/environment.txt")
reticulate::py_install(packages = required_pkgs, envname = env_name, method = "virtualenv")

# Confirm active Python path
py_config()$python
```

#### Example usage:

```r
run_PRIMEloci_example(python_path = py_config()$python)
```

---

### 🌬️ Option 2: Conda

`environment.yml` is included in the `inst/envfile` folder.

#### How to create the environment:

```bash
cd PRIME

# Recommended: NumPy 2.x
conda env create -f inst/envfile/environment.yml

# Activate the environment
conda activate prime-env
which python3
conda deactivate

# Copy this path to use as python_path in R when calling run_PRIMEloci_example() or run_PRIMEloci_facet_example()
```

#### Example in R:

```r
run_PRIMEloci_example(python_path = "~/.conda/envs/prime-env/bin/python3")
run_PRIMEloci_facet_example(python_path = "~/.conda/envs/prime-env/bin/python3")
```

---

### 🧪 Option 3: Virtualenv (manual setup)

This is an advanced manual option for setting up the environment outside of R. It is useful if managing environments via shell or external tools.

#### Setup instructions

```bash
python3 -m venv ~/prime_env
source ~/prime_env/bin/activate
pip3 install --upgrade pip
pip3 install -r inst/envfile/environment.txt
which python3

# Copy this path to use as python_path in R when calling run_PRIMEloci_example() or run_PRIMEloci_facet_example()
```

#### Example in R:

```r
run_PRIMEloci_example(python_path = "~/prime_env/bin/python3")
run_PRIMEloci_facet_example(python_path = "~/prime_env/bin/python3")
```

---

### 🔧 Option 4: Use existing Python installation with `requirements.txt`

If you already have a working Python installation and want to use it directly:

```bash
pip3 install --upgrade pip
pip3 install -r inst/envfile/environment.txt
which python3

# Use this path in your R function call
```

```r
run_PRIMEloci_example(python_path = "/path/to/your/python")
```

---

## 🧠 Tips

- Always check the current Python path with `which python3` **after** activating your environment.
- Use that full path in the `python_path` argument in R.
- All four environments (reticulate virtualenv, conda, manual virtualenv, and existing Python) are compatible with `PRIMEloci()`.

---

© 2025 PRIME setup protocol | Supports PRIMEloci-compatible Python integration
