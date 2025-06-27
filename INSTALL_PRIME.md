# 🧪 PRIME Installation Preparation Guide

This guide explains how to set up the required **R** and **Python** environments to use the `PRIME` R package. The Python environment is also compatible with `PRIMEloci()` if used.

```bash
# In terminal
git clone ⁦https://github.com/anderssonlab/PRIME.git
```

FOR MacOS users: **`libomp` is required** for LightGBM to enable OpenMP (multithreading).

```bash
# Check if libomp exists

## Apple Silicon (M1/M2/M3/...)
ls /opt/homebrew/opt/libomp/lib/libomp.dylib
## Intel macOS (x86_64)
ls /usr/local/opt/libomp/lib/libomp.dylib
```

```bash
# If libomp does not exist

brew install libomp

## If Apple Silicon (M1/M2/M3/...)
echo 'export DYLD_LIBRARY_PATH="/opt/homebrew/opt/libomp/lib:$DYLD_LIBRARY_PATH"' >> ~/.zshrc
## or if Intel macOS (x86_64)
echo 'export DYLD_LIBRARY_PATH="/usr/local/opt/libomp/lib:$DYLD_LIBRARY_PATH"' >> ~/.zshrc

source ~/.zshrc
```

---

## 📦 R Installation for PRIME

The R package `PRIME` depends on a mix of CRAN and Bioconductor packages, and a custom GitHub version of `bcp`.

### ✅ Full R Setup

```r
# 1. Install required CRAN packages
install.packages(c(
	"R.utils",
  "future",
  "future.apply",
  "future.callr",
  "foreach",
  "argparse",
  "doParallel",
  "reticulate",
  "arrow",
  "Matrix",
  "caTools",
  "zoo"
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

# If installing bcp on macOS fails due to missing gcc,
# check with << gfortran --version >>. 
# If not found, install Xcode command line tools and Homebrew.
# Then run << brew install gcc >>, check << gfortran --version >>, 
# and reinstall bcp.

# Using remotes package with << remote::install_github("swang87/bcp") >>
# instead of devtools if not successfully install devtools

# 7. Install PRIME
## Install from .tar.gz (model will be set up at PRIME inst directory),
## otherwise, make sure that the model in the path was exist.
install.packages("/PATH/TO/PRIME/PRIME_0.1.1.6.tar.gz")
```

---

## 🐍 Python Environment Setup (Compatible with PRIMEloci functions)

You can prepare the Python environment in **four ways**. The recommended default is the `reticulate`-managed virtualenv using **NumPy 2.x** (use NumPy 1.x only if required for legacy compatibility)

---

### 🌐 Option 1 : Virtualenv via `reticulate` in R

This is the **default and recommended** method for setting up the Python environment using only R. It ensures full compatibility with `PRIME` and `PRIMEloci()` 

Reticulate virtualenv (in R) is easiest for R-focused workflows without needing external software, **but requires running `use_virtualenv()` in each session** and is not ideal for CLI use.

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
required_pkgs <- readLines("/PATH/TO/PRIME/inst/envfile/environment.txt")
reticulate::py_install(packages = required_pkgs, envname = env_name, method = "virtualenv")

# Confirm active Python path
py_config()$python
```

#### Example usage:

```r
run_PRIMEloci_focal_example(python_path = py_config()$python)
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
run_PRIMEloci_focal_example(python_path = "~/.conda/envs/prime-env/bin/python3")
run_PRIMEloci_example(python_path = "~/.conda/envs/prime-env/bin/python3")
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
run_PRIMEloci_focal_example(python_path = "~/prime_env/bin/python3")
run_PRIMEloci_example(python_path = "~/prime_env/bin/python3")
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
run_PRIMEloci_focal_example(python_path = "/path/to/your/python")
run_PRIMEloci_example(python_path = "/path/to/your/python")
```

---

## 🧠 Tips

- Always check the current Python path with `which python3` **after** activating your environment.
- Use that full path in the `python_path` argument in R.
- All four environments (reticulate virtualenv, conda, manual virtualenv, and existing Python) are compatible with `PRIMEloci()`.

---

© 2025 PRIME setup protocol | Supports PRIMEloci-compatible Python integration
