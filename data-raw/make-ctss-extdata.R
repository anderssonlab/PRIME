# Regenerate packaged CTSS example objects for vignettes.
# Run from the package root with:
#   Rscript data-raw/make-ctss-extdata.R

suppressPackageStartupMessages({
  library(CAGEfightR)
  library(GenomeInfoDb)
  library(PRIME)
  library(rtracklayer)
  library(SummarizedExperiment)
})

find_package_root <- function(starts) {
  starts <- unique(normalizePath(starts, mustWork = FALSE))
  starts <- starts[!is.na(starts) & nzchar(starts)]

  for (start in starts) {
    path <- if (dir.exists(start)) start else dirname(start)
    repeat {
      if (file.exists(file.path(path, "DESCRIPTION"))) {
        return(path)
      }
      parent <- dirname(path)
      if (identical(parent, path)) {
        break
      }
      path <- parent
    }
  }

  stop("Could not find package root; run from within the PRIME checkout.")
}

with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(expr)
}

file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
frame_files <- vapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) NA_character_ else frame$ofile
}, character(1))
package_root <- find_package_root(c(file_arg, frame_files, getwd()))

design_file <- file.path(package_root, "inst", "extdata", "design_matrix_first10pct.tsv")
bw_dir <- file.path(package_root, "inst", "extdata", "cage_bw")
out_ctss <- file.path(package_root, "inst", "extdata", "ctss.rds")
out_ctss_clean <- file.path(package_root, "inst", "extdata", "ctss_clean.rds")

if (!file.exists(design_file)) {
  stop("Missing design matrix: ", design_file)
}
if (!dir.exists(bw_dir)) {
  stop("Missing BigWig directory: ", bw_dir)
}

message("Reading design matrix: ", design_file)
design <- read.table(design_file, header = TRUE, sep = "\t")
rownames(design) <- design$Name

bw_plus <- rtracklayer::BigWigFileList(file.path(bw_dir, design$BigWigPlus))
bw_minus <- rtracklayer::BigWigFileList(file.path(bw_dir, design$BigWigMinus))
names(bw_plus) <- names(bw_minus) <- design$Name

message("Quantifying CTSSs")
ctss <- CAGEfightR::quantifyCTSSs(
  plusStrand = bw_plus,
  minusStrand = bw_minus,
  design = design
)
ctss <- CAGEfightR::calcTotalTags(ctss)
ctss <- CAGEfightR::calcTPM(ctss)
ctss <- CAGEfightR::calcPooled(ctss)
ctss <- GenomeInfoDb::keepStandardChromosomes(ctss, pruning.mode = "coarse")

message("Creating singleton-filtered CTSS object")
min_depth <- min(SummarizedExperiment::colData(ctss)$totalTags)
# Keep regenerated RDS files reproducible without leaking RNG changes to
# interactive sessions that source this script.
ctss_sub <- with_seed(1, PRIME::subsampleTarget(ctss, target = min_depth))
ctss_sub <- CAGEfightR::calcTPM(ctss_sub)
ctss_clean <- PRIME::rmSingletons(ctss_sub)

message("Writing: ", out_ctss)
saveRDS(ctss, out_ctss, compress = "xz")
message("Writing: ", out_ctss_clean)
saveRDS(ctss_clean, out_ctss_clean, compress = "xz")

message("Done")
