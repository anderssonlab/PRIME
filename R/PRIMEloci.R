# profile_log_message("Step 3 completed", "my_pipeline.log")
# profile_log_message("⚠️ Skipping chromosome chrX", "my_pipeline.log", level = "WARN")
# profile_log_message("❌ Failed to load file", "my_pipeline.log", level = "ERROR")
# log_file <- "my_pipeline.log"
# ✅
# profile_log_message("🚀 Pipeline started", log_file)
# profile_log_message("🔍 Processing chr1", log_file)
# 
# 
# log_file <- file.path(log_dir, "PRIMEloci-2.log")
# 
# plc_log("Started TC validation", log_file)
# plc_log(sprintf("⚠️ Skipping chromosome %s", chr_name), log_file, level = "WARN")
# plc_log("TC validation complete", log_file, print_console = FALSE)

#' PRIMEloci function
#'
#' @export
#'
#' @import GenomicRanges
#' @import assertthat
PRIMEloci <- function(ctss_rse,
                      tc_object = NULL,
                      outdir,
                      ext_dis = 200,
                      save_tc = TRUE,
                      tc_object_name = "tc_grl.rds",
                      sld_window = TRUE,
                      sld_by = 20,
                      sld_object_name = "sld_tc_grl.rds",
                      keep_tmp = TRUE,
                      num_cores = NULL) {

  message("Running PRIMEloci ...")

  # Ensure directories exist
  if (dir.exists(outdir)) {
    warning("Warning: Output directory '", outdir,
            "' already exists. !!! Files may be overwritten. !!!")
  } else {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  primeloci_tmp <- file.path(outdir, "PRIMEloci_tmp")
  dir.create(primeloci_tmp, recursive = TRUE, showWarnings = FALSE)

  log_dir <- file.path(outdir, "PRIMEloci_log")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  # Ensure numeric parameters are integers
  assertthat::assert_that(is.numeric(ext_dis),
                          msg = "ext_dis must be numeric.")
  assertthat::assert_that(ext_dis %% 1 == 0,
                          msg = "ext_dis must be an integer.")

  assertthat::assert_that(is.numeric(sld_by), msg = "sld_by must be numeric.")
  assertthat::assert_that(sld_by %% 1 == 0, msg = "sld_by must be an integer.")

  # _2_
  log_2 <- file.path(log_dir, "PRIMEloci-2.log")
  plc_log("\n\n\n 🚀 Running PRIMEloci -2 : get extended tc or validate the tc object provided", log_2) # nolint: line_length_linter.

  if (is.null(tc_object)) {
    plc_log("🔍 Creating tc object\n", log_2)
    tc_grl <- plc_trylog(
      get_tcs_and_extend_fromthick(ctss_rse, ext_dis = ext_dis),
      log_file = log_2,
      print_console = FALSE
    )
    if (save_tc) saveRDS(tc_grl, file = file.path(outdir, tc_object_name))
  } else {
    tc_grl <- tc_object
  }

  plc_log("🔍 Validating tc object\n", log_2)
  validate_tc <- plc_trylog(
    validate_tc_object(tc_grl, ctss_rse, ext_dis = ext_dis),
    log_file = log_2,
    print_console = TRUE
  )

  if (!validate_tc) {
    msg <- "\n❌ TC object validation failed. Ensure the TC object is valid."
    plc_log(msg, log_2, level = "ERROR")
    stop(msg)
  }
  plc_log("✅ DONE :: TC object is validated and ready to use.", log_2)

  # _3_
  if (sld_window) {
    log_3 <- file.path(log_dir, "PRIMEloci-3.log")
    plc_log("\n\n\n 🚀 Running PRIMEloci -3 : sliding windows covering reduced TC regions", log_3) # nolint: line_length_linter.

    if (inherits(tc_grl, "GenomicRanges::GRanges")) {
      start_time <- Sys.time()
      plc_log("Processing single GRanges object", log_3)
      tc_sliding_window_grl <- plc_trylog({
        tc_sliding_window(tc_grl,
                          sld_by = sld_by,
                          ext_dis = ext_dis,
                          num_cores = num_cores)
      }, log_file = log_3)
      plc_log(sprintf("Time taken: %.2f minutes",
                      as.numeric(difftime(Sys.time(),
                                          start_time,
                                          units = "mins"))),
              log_3)
    } else if (inherits(tc_grl, "GenomicRanges::GRangesList") ||
                 inherits(tc_grl, "CompressedGRangesList")) {
      tc_sliding_window_grl <- lapply(seq_along(tc_grl), function(i) {
        start_time <- Sys.time()
        gr_name <- if (!is.null(names(tc_grl))) names(tc_grl)[i] else paste0("Sample_", i) # nolint: line_length_linter.
        plc_log(sprintf("Processing: %s", gr_name), log_3)
        result <- tc_sliding_window(tc_grl[[i]],
                                    sld_by = sld_by,
                                    ext_dis = ext_dis,
                                    log_file = log_3,
                                    num_cores = num_cores)
        plc_log(sprintf("Time taken: %.2f minutes",
                        as.numeric(difftime(Sys.time(),
                                            start_time,
                                            units = "mins"))),
                log_3)
        return(result)
      })
      if (!is.null(tc_sliding_window_grl) &&
            length(tc_sliding_window_grl) > 0) {
        tc_sliding_window_grl <- GenomicRanges::GRangesList(tc_sliding_window_grl) # nolint: line_length_linter.
      } else {
        msg <- "❌ Processed TC object list is empty. Ensure tc_grl contains valid data." # nolint: line_length_linter.
        plc_log(msg, log_3, level = "ERROR")
        stop(msg)
      }

    } else {
      msg <- "❌ tc_grl must be either a GRanges, GRangesList, or CompressedGRangesList object." # nolint: line_length_linter.
      plc_log(msg, log_3, level = "ERROR")
      stop(msg)
    }

    plc_log(sprintf("Saving TC objects to PRIMEloci_tmp .."),
            log_3)
    saveRDS(tc_sliding_window_grl, file.path(primeloci_tmp, sld_object_name))
    plc_log("✅ DONE :: Sliding window TC object is saved to PRIMEloci_tmp.",
            log_3)
    tc_for_profile <- tc_sliding_window_grl
  } else {
    tc_for_profile <- tc_grl
  }

  assertthat::assert_that(!is.null(tc_for_profile),
                          msg = "❌ tc_for_profile is NULL at return. Something failed.") # nolint: line_length_linter.
  return(tc_for_profile)

}
