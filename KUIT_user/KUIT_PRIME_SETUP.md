# ⚙️ PRIME Installation on KUIT Server

This guide outlines a full PRIME and PRIMEloci-compatible environment setup on a KUIT server.

---

## 🧩 Step 1: Load Required Modules

```bash
module load gcc/13.2.0 openjdk/20.0.0 R/4.4.0
module load hdf5 netcdf-c/4.8.1
module load miniconda/24.5.0
```

---

## 📥 Step 2: Clone PRIME Repository

```bash
git clone https://github.com/anderssonlab/PRIME.git
cd PRIME
```

---

## 📦 Step 3: Install PRIME R Package

```bash
R

# install necessary packages 
# see /PRIME/PREP_INSTALL_PRIME.md

install.packages("/PATH/TO/PRIME/PRIME_0.1.1.6.tar.gz", lib="/PATH/TO/R_packakges/")
q()
```

---

## 🐍 Step 4: Python Conda Environment Setup

```bash
# Check Python path if already installed
which python3

# Suggested method: Use local conda env (NumPy 2.x preferred)
cd PRIME/inst/envfile
conda env create -f environment.yml
conda activate prime-env 
which python3

# Output example:
# ~/.conda/envs/prime-env/bin/python3
```

---

## 🔗 Step 5: Export libcurl Settings (Optional but Recommended)

If you encounter libcurl-related issues when using R with conda Python:

- With in conda env

```bash

export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export LD_PRELOAD=$CONDA_PREFIX/lib/libcurl.so
R
```

- Outside conda env

```bash
# Check
ls $HOME/.conda/envs/prime-env/lib/libcurl*

# Expected output includes:
# `libcurl.so` (symlink)
# `libcurl.so.4` (symlink)
# `libcurl.so.4.8.0` (actual binary)
```

```bash
export LD_LIBRARY_PATH=/home/KU_ID/.conda/envs/prime-env/lib:$LD_LIBRARY_PATH
export LD_PRELOAD=/home/KU_ID/.conda/envs/prime-env/lib/libcurl.so
R
```

---

## 🚀 Step 6: Start running in R

```R
# Set .libPaths if necessary
.libPaths("/path/to/R_packakges/")

# Set env to ignore ~/.local packages, if conda was set
Sys.setenv(PYTHONNOUSERSITE = "1")

library(GenomicRanges)
library(PRIME)

# Run examples
plc_facet <- run_PRIMEloci_focal_example("~/.conda/envs/prime-env/bin/python3")
plc <- run_PRIMEloci_example("~/.conda/envs/prime-env/bin/python3")
```

---

© 2025 PRIME setup protocol for KUIT server
