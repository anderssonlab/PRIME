library(SummarizedExperiment)
library(GenomicRanges)
library(PRIME)

# Load your real CTSS RangedSummarizedExperiment
ctss_rse <- readRDS("/Users/natsudanav/Desktop/PRIME/tests/testdata/ctss_rse.rds") # nolint

# Output directory for this test
test_outdir <- "/Users/natsudanav/Desktop/PRIME/tests/testoutput"
dir.create(test_outdir, showWarnings = FALSE, recursive = TRUE)

PRIME::PRIMEloci(
  ctss_rse,
  tc_object = NULL,
  test_outdir,
  save_tc = TRUE,
  tc_object_name = "tc_grl.rds",
  sld_window = TRUE,
  sld_object_name = "sld_tc_grl.rds",
  sld_by = 20,
  profile_dir_name = "PRIMEloci_profiles",
  file_type = "npz",
  addtn_to_filename = "test_",
  save_count_profiles = TRUE,
  python_path = "~/.virtualenvs/prime-env",
  keep_tmp = TRUE,
  num_cores = NULL,
  ext_dis = 200
)

cat("✅ Test with real ctss_rse complete.\nCheck output in:", test_outdir, "\n")
