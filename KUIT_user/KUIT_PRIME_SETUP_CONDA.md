# PRIME Environment Setup Guide (with Expected Outputs)

This guide sets up PRIME in a clean, reproducible `conda` environment.

---

## 1. Clean start
```bash
module purge
module load miniconda/24.5.0
conda config --set channel_priority strict
```

**Expect:** No errors. `conda config` confirms channel priority is set.

---

## 2. Remove old env (if exists) and create new one
```bash
conda env remove -n prime-conda-env -y >/dev/null 2>&1 || true
conda create -y -n prime-conda-env -c conda-forge python=3.11 r-base=4.4.*
```

**Expect:** Conda solves environment and installs Python 3.11.x + R 4.4.x.

---

## 3. Activate & sanitize
```bash
conda activate prime-conda-env
unset LD_PRELOAD
unset LD_LIBRARY_PATH
```

**Expect:** Prompt changes to `(prime-conda-env)`.

---

## 4. Ignore `~/.local` Python site-packages inside this env
```bash
conda env config vars set PYTHONNOUSERSITE=1
conda deactivate && conda activate prime-conda-env
```

**Expect:** Env reactivates cleanly. No `.local` packages will be seen.

---

## 5. Make R prefer this env’s library automatically
```bash
RPROFILE_SITE="$CONDA_PREFIX/lib/R/etc/Rprofile.site"
mkdir -p "$(dirname "$RPROFILE_SITE")"
cat > "$RPROFILE_SITE" <<'RPROF'
## Prefer the conda env's library first
conda_lib <- "~/.conda/envs/prime-conda-env/lib/R/library"
if (dir.exists(path.expand(conda_lib))) {
  .libPaths(unique(c(path.expand(conda_lib), .libPaths())))
}
## Avoid user lib overriding this env
Sys.unsetenv("R_LIBS_USER")
RPROF
```

Verify:
```bash
R -q -e "print(.libPaths())"
```

**Expect:** First entry is `~/.conda/envs/prime-conda-env/lib/R/library`.

---

## 6. Install Python dependencies
```bash
## Change path to PRIME directiry before running this command
conda install -y -c conda-forge --update-specs --file /PATH/TO/PRIME/inst/envfile/environment.txt
```

**Expect:** Conda installs specified versions.

Quick check:
```bash
python - <<'PY'
import sys, numpy, scipy, pandas, sklearn, joblib, pyarrow, fastparquet, lightgbm
print("python", sys.version.split()[0])     # expect 3.11.x
print("numpy", numpy.__version__)           # expect 2.0.2
print("scipy", scipy.__version__)           # expect 1.14.x
print("pandas", pandas.__version__)         # expect 2.2.x (within >=2.1,<2.3)
print("scikit-learn", sklearn.__version__)  # expect 1.6.x (within >=1.4,<1.7)
print("joblib", joblib.__version__)         # expect >=1.2
print("pyarrow", pyarrow.__version__)       # expect <18
print("fastparquet", fastparquet.__version__) # expect >=2023.4
print("lightgbm", lightgbm.__version__)     # expect 4.6.0
PY
```

---

## 7. Install prebuilt CRAN R packages
```bash
conda install -y -c conda-forge   r-r.utils r-future r-future.apply r-future.callr r-foreach r-argparse   r-doparallel r-reticulate r-arrow r-igraph r-catools r-zoo   r-biocmanager r-remotes r-devtools
```

**Expect:** R packages installed without compilation.

---

## 8. Install Bioconductor packages (Conda binaries)
```bash
conda install -y -c conda-forge -c bioconda   bioconductor-cagefightr   bioconductor-rtracklayer   bioconductor-genomicranges   bioconductor-iranges   bioconductor-genomeinfodb   bioconductor-summarizedexperiment   bioconductor-biocparallel   bioconductor-bsgenome
```

Quick check:
```bash
R -q -e 'library(CAGEfightR); packageVersion("CAGEfightR")'
```

**Expect:** Prints a version (e.g., ‘1.26.0’) with no segfaults.

---

## 9. Install compilers (for building PRIME)
```bash
conda install -y -c conda-forge gxx_linux-64 gcc_linux-64 gfortran_linux-64 make pkg-config
```

Verify inside R:
```r
R -q -e 'Sys.getenv(c("CC","CXX","FC")); system("R CMD config CC"); system("R CMD config CXX"); system("R CMD config FC")'
```

**Expect:** Compiler paths inside conda, like:
```
x86_64-conda-linux-gnu-cc
x86_64-conda-linux-gnu-c++ -std=gnu++17
x86_64-conda-linux-gnu-gfortran
```

---

## 10. (Optional) Install bcp from GitHub
```bash
R -q -e 'Sys.unsetenv("R_LIBS_USER"); target_lib <- .libPaths()[1];
if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", lib=target_lib, repos="https://cloud.r-project.org");
remotes::install_github("swang87/bcp", lib=target_lib, upgrade="never", dependencies=TRUE, force=TRUE);
library(bcp, lib.loc=target_lib); cat("bcp loaded from: ", system.file(package="bcp"), "\n"); print(packageVersion("bcp"))'
```

**Expect:** `bcp` loads and version prints.

---

## 11. Install PRIME (local tarball)
```bash
## Change path to PRIME directiry before running this command
R -q -e 'Sys.setenv(RETICULATE_PYTHON="~/.conda/envs/prime-conda-env/bin/python3", PYTHONNOUSERSITE="1"); install.packages("/PATH/TO/PRIME/PRIME_0.1.1.6.tar.gz", repos=NULL, type="source", lib=.libPaths()[1]); library(PRIME); packageVersion("PRIME"); library(reticulate); print(py_config()); cat("\nPRIME loaded OK ✅\n")'
```

**Expect:**
- PRIME installs cleanly.
- `packageVersion("PRIME")` prints `0.1.1.6`.
- `py_config()` shows Python from `prime-conda-env`.
- Final message: `PRIME loaded OK ✅`.

---

## 12. Quick PRIME test
```r
R
library(GenomicRanges)
library(PRIME)
packageVersion("PRIME")
plc_focal_example <- run_PRIMEloci_focal_example(python_path = "~/.conda/envs/prime-conda-env/bin/python3")
plc_example <- run_PRIMEloci_example(python_path = "~/.conda/envs/prime-conda-env/bin/python3")
```

**Expect:** Example runs without error.

---

✅ Done — Environment `prime-conda-env` is ready.
