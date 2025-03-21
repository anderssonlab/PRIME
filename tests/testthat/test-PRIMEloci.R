test_that("PRIMEloci returns correct object type", {
  ctss_rse <- readRDS("/Users/natsudanav/Desktop/PRIME/tests/testdata/ctss_rse.rds") # nolint: line_length_linter.

  # Case 1: Sliding window is TRUE (expect GRangesList)
  result1 <- PRIMEloci(ctss_rse, outdir = tempdir(), sld_window = TRUE)
  expect_true(inherits(result1, "GRangesList") ||
                inherits(result1, "CompressedGRangesList"))

  # Case 2: Sliding window is FALSE (expect same type as Step 2)
  result2 <- PRIMEloci(ctss_rse, outdir = tempdir(), sld_window = FALSE)

  # Expect same type as tc_grl in Step 2 (allow GRanges or GRangesList)
  expect_true(inherits(result2, "GenomicRanges::GRanges") ||
                inherits(result2, "GenomicRanges::GRangesList") ||
                inherits(result2, "CompressedGRangesList"))
})
