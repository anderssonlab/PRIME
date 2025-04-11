# 🧪 PRIME Installation Preparation Guide

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
    "FC = /usr/local/bin/gfortran
",
    "F77 = /usr/local/bin/gfortran
",
    "FLIBS = -L/usr/local/lib/gcc/13 -lgfortran -lquadmath
",
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
), lib = "/path/to/R_packakges/")

# 5. Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools", lib = "/path/to/R_packakges/"")

# 6. Install bcp from GitHub (custom version)
devtools::install_github("swang87/bcp")

# ****** remotes instead of devtools if not successfully install devtools
```

---

## 🐍 Python Environment Setup (Compatible with PRIMEloci)

You can prepare the Python environment in **four ways**. The recommended default is the `reticulate`-managed virtualenv using **NumPy 2.x**.

---

### 🌐 Option 1 (Default): Virtualenv via `reticulate` in R (NumPy 2.x)

This is the **default and recommended** method for setting up the Python environment using only R. It ensures full compatibility with `PRIME` and `PRIMEloci()` using NumPy 2.x.

#### Required Python packages (NumPy 2.x version):

```r
required_pkgs <- c(
  "numpy==2.0.0",
  "scipy==1.11.4",
  "pandas==2.1.4",
  "joblib==1.3.2",
  "scikit-learn==1.3.2",
  "lightgbm==4.1.0",
  "pyarrow==14.0.2",
  "fastparquet==2023.8.0"
)
```

#### Setup instructions in R:

```r
library(reticulate)

# Define the environment name
env_name <- "prime_env_np2"

# Create the virtual environment if it doesn't exist
virtualenv_create(envname = env_name)

# Use and configure it
use_virtualenv(env_name, required = TRUE)

# Install Python packages
py_install(packages = required_pkgs, envname = env_name, method = "virtualenv")

# Confirm active Python path
py_config()$python
```

#### Example usage:

```r
run_PRIMEloci_example(python_path = py_config()$python)
```

---

### 🌬️ Option 2: Conda (alternative setup with optional NumPy 1.x or 2.x)

> 💡 **Note:** NumPy 2.x is recommended. Use NumPy 1.x only if required for legacy compatibility.

`environment_np2.yml` (recommended) and `environment_np1.yml` (optional legacy) are included in the `inst/` folder.

#### How to create the environment:

```bash
cd PRIME

# Recommended: NumPy 2.x
dconda env create -f inst/conda_env/environment_np2.yml

# Optional legacy: NumPy 1.x (only if needed)
# conda env create -f inst/conda_env/environment_np1.yml

# Activate the environment
conda activate prime-env-numpy2
which python
# Copy this path to use as python_path in R when calling run_PRIMEloci_example() or run_PRIMEloci_facet_example()
```

#### Example in R:

```r
run_PRIMEloci_example(python_path = "~/.conda/envs/prime-env-numpy2/bin/python")
run_PRIMEloci_facet_example(python_path = "~/.conda/envs/prime-env-numpy2/bin/python")
```

---

### 🧪 Option 3: Virtualenv (manual setup)

This is an advanced manual option for setting up the environment outside of R. It is useful if managing environments via shell or external tools.

#### `requirements.txt`

```txt
numpy==2.0.0
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
python3 -m venv ~/prime_env_np2
source ~/prime_env_np2/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
which python
# Copy this path to use as python_path in R when calling run_PRIMEloci_example() or run_PRIMEloci_facet_example()
```

#### Example in R:

```r
run_PRIMEloci_example(python_path = "~/prime_env_np2/bin/python")
run_PRIMEloci_facet_example(python_path = "~/prime_env_np2/bin/python")
```

---

### 🔧 Option 4: Use existing Python installation with `requirements.txt`

If you already have a working Python installation and want to use it directly:

```bash
pip install --upgrade pip
pip install -r requirements.txt
which python
# Use this path in your R function call
```

```r
run_PRIMEloci_example(python_path = "/path/to/your/python")
```

---

## 🧠 Tips

- Always check the current Python path with `which python` **after** activating your environment.
- Use that full path in the `python_path` argument in R.
- All four environments (reticulate virtualenv, conda, manual virtualenv, and existing Python) are compatible with `PRIMEloci()`.

---

© 2025 PRIME setup protocol | Supports PRIMEloci-compatible Python integration
