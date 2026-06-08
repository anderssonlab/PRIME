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

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (is.na(script_file)) {
  script_file <- file.path("data-raw", "make-ctss-extdata.R")
}
package_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
if (!file.exists(file.path(package_root, "DESCRIPTION"))) {
  package_root <- normalizePath(getwd(), mustWork = TRUE)
}

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
set.seed(1)
min_depth <- min(SummarizedExperiment::colData(ctss)$totalTags)
ctss_sub <- PRIME::subsampleTarget(ctss, target = min_depth)
ctss_sub <- CAGEfightR::calcTPM(ctss_sub)
ctss_clean <- PRIME::rmSingletons(ctss_sub)

message("Writing: ", out_ctss)
saveRDS(ctss, out_ctss, compress = "xz")
message("Writing: ", out_ctss_clean)
saveRDS(ctss_clean, out_ctss_clean, compress = "xz")

message("Done")
