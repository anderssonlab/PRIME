#' PRIMEloci function
#'
#' @export
#'
#' @import GenomicRanges
#' @import assertthat
PRIMEloci <- function(
    ctss_rse,
    outdir,
    python_path = "~/.virtualenvs/prime-env",
    score_threshold = 0.75,
    score_diff = 0.1,
    num_cores = NULL,
    keep_tmp = FALSE) {

  # Validate inputs

  # Check ctss_rse
  assertthat::assert_that(
    methods::is(ctss_rse, "RangedSummarizedExperiment"),
    msg = "`ctss_rse` must be a RangedSummarizedExperiment object."
  )

  # Normalize and verify output directory
  outdir <- normalizePath(outdir, mustWork = FALSE)
  if (dir.exists(outdir)) {
    warning("⚠️ Output directory already exists: ", outdir)
  } else {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  outdir <- normalizePath(outdir, mustWork = TRUE)

  log_dir <- file.path(outdir, "PRIMEloci_log")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  primeloci_tmp <- file.path(outdir, "PRIMEloci_tmp")
  dir.create(primeloci_tmp, recursive = TRUE, showWarnings = FALSE)

  #Check numeric parameters
  assertthat::assert_that(
    is.numeric(score_threshold),
    score_threshold > 0,
    score_threshold < 1,
    msg = "`score_threshold` must be a numeric value between 0 and 1."
  )

  assertthat::assert_that(
    is.numeric(score_diff),
    score_diff >= 0,
    score_diff < score_threshold,
    msg = "`score_diff` must be a non-negative numeric value and smaller than `score_threshold`." # nolint: line_length_linter.
  )

  if (!is.null(num_cores)) {
    assertthat::assert_that(
      is.numeric(num_cores),
      num_cores %% 1 == 0,
      num_cores > 0,
      msg = "`num_cores` must be a positive integer or NULL."
    )
  }

  log_1 <- file.path(log_dir, "PRIMEloci-1.log")

  plc_trylog({

    plc_log("🚀 Starting PRIMEloci pipeline", log_1)

    # Set up Python environment
    python_path <- path.expand(python_path)
    reticulate::use_virtualenv(python_path, required = TRUE)
    py_conf <- reticulate::py_config()

    # keep version info
    plc_log(sprintf("🔧 R version: %s", R.version.string),
            log_1, print_console = FALSE)
    plc_log(sprintf("📦 Loaded Python: %s", py_conf$python),
            log_1, print_console = FALSE)
    lapply(capture.output(py_conf), function(line) plc_log(paste("📦", line), log_1, print_console = FALSE)) # nolint: line_length_linter.

  }, log_file = log_1)

  # setting
  tc_object <- NULL
  save_tc <- FALSE
  tc_object_name <- "tc_grl.rds"
  sld_by <- 20
  sld_object_name <- "sld_tc_grl.rds"
  save_sld_tc <- FALSE
  profile_dir_name <- "PRIMEloci_profiles"
  postprocess_partial_name <- "pred_all"
  save_count_profiles <- FALSE
  ext_dis <- 200
  file_type <- "npz"
  addtn_to_filename <- ""
  name_prefix <- "PRIMEloci"
  model_name <- "PRIMEloci_GM12878_model_1.0.sav"
  core_width <- 151


  # _2_
  log_2 <- file.path(log_dir, "PRIMEloci-2.log")
  plc_log("\n\n\n 🚀 Running PRIMEloci -2 : get extended tc or validate the tc object provided", log_2) # nolint: line_length_linter.

  if (is.null(tc_object)) {
    plc_log("🔹 Creating tc object\n", log_2)
    tc_grl <- plc_trylog(
      get_tcs_and_extend_fromthick(ctss_rse, ext_dis = ext_dis),
      log_file = log_2,
      print_console = FALSE
    )
    if (save_tc) saveRDS(tc_grl, file = file.path(outdir, tc_object_name))
  } else {
    tc_grl <- tc_object
  }

  plc_log("🔹 Validating tc object\n", log_2)
  validate_tc <- plc_trylog(
    validate_tc_object(tc_grl, ctss_rse, ext_dis = ext_dis),
    log_file = log_2,
    print_console = TRUE
  )

  if (!validate_tc) {
    msg <- "\nTC object validation failed. Ensure the TC object is valid."
    plc_log(msg, log_2, level = "❌ ERROR")
    stop(msg)
  }
  plc_log("✅ DONE :: TC object is validated and ready to use.", log_2)


  # _3_
  log_3 <- file.path(log_dir, "PRIMEloci-3.log")
  plc_log("\n\n\n 🚀 Running PRIMEloci -3 : sliding windows covering reduced TC regions", log_3) # nolint: line_length_linter.

  if (inherits(tc_grl, "GenomicRanges::GRanges")) {
    start_time <- Sys.time()
    plc_log("🔹 Processing single GRanges object", log_3)
    tc_sliding_window_grl <- plc_trylog({
      tc_sliding_window(tc_grl,
                        sld_by = sld_by,
                        ext_dis = ext_dis,
                        num_cores = num_cores)
    }, log_file = log_3)
    plc_log(sprintf("⏱️ Time taken: %.2f minutes",
                    as.numeric(difftime(Sys.time(),
                                        start_time,
                                        units = "mins"))),
            log_3)
  } else if (inherits(tc_grl, "GenomicRanges::GRangesList") ||
               inherits(tc_grl, "CompressedGRangesList")) {
    tc_sliding_window_grl <- lapply(seq_along(tc_grl), function(i) {
      start_time <- Sys.time()
      gr_name <- if (!is.null(names(tc_grl))) names(tc_grl)[i] else paste0("Sample_", i) # nolint: line_length_linter.
      plc_log(sprintf("🔹 Processing: %s", gr_name), log_3)
      result <- tc_sliding_window(tc_grl[[i]],
                                  sld_by = sld_by,
                                  ext_dis = ext_dis,
                                  log_file = log_3,
                                  num_cores = num_cores)
      plc_log(sprintf("⏱️ Time taken: %.2f minutes",
                      as.numeric(difftime(Sys.time(),
                                          start_time,
                                          units = "mins"))),
              log_3)
      result
    })

    if (!is.null(tc_sliding_window_grl) &&
          length(tc_sliding_window_grl) > 0) {
      tc_sliding_window_grl <- GenomicRanges::GRangesList(tc_sliding_window_grl) # nolint: line_length_linter.
      names(tc_sliding_window_grl) <- names(tc_grl)
    } else {
      msg <- "Processed TC object list is empty. Ensure tc_grl contains valid data." # nolint: line_length_linter.
      plc_log(msg, log_3, level = "❌ ERROR")
      stop(msg)
    }
  } else {
    msg <- "tc_grl must be either a GRanges, GRangesList, or CompressedGRangesList object." # nolint: line_length_linter.
    plc_log(msg, log_3, level = "❌ ERROR")
    stop(msg)
  }

  if (save_sld_tc) {
    plc_log(sprintf("Saving TC objects to PRIMEloci_tmp .."),
            log_3)
    saveRDS(tc_sliding_window_grl,
            file.path(primeloci_tmp, sld_object_name))
  }
  plc_log("✅ DONE :: Sliding window TC object is saved to PRIMEloci_tmp.",
          log_3)
  tc_for_profile <- tc_sliding_window_grl

  assertthat::assert_that(!is.null(tc_for_profile),
                          msg = "❌ tc_for_profile is NULL at return. Something failed.") # nolint: line_length_linter.


  # _4_
  log_file_name_4 <- "PRIMEloci-4.log"
  log_4 <- file.path(log_dir, log_file_name_4)
  plc_log("\n\n\n 🚀 Running PRIMEloci -4 : compute count & normalized profiles for each sample", # nolint: line_length_linter.
          log_4)

  plc_trylog(
    PRIMEloci_profile(
      ctss_rse,
      tc_for_profile,
      outdir,
      profile_dir_name,
      file_type = file_type,
      addtn_to_filename = addtn_to_filename,
      save_count_profiles = save_count_profiles,
      num_cores = num_cores,
      log_file = log_4,
      ext_dis
    ),
    log_4
  )


  # _5_
  log_5 <- file.path(log_dir, "PRIMEloci-5.log")
  plc_log("\n\n\n 🚀 Running PRIMEloci -5 : Prediction using PRIMEloci model",
          log_5)

  profile_main_dir <- file.path(primeloci_tmp,
                                profile_dir_name)

  profiles_subtnorm_dir <- file.path(profile_main_dir, "profiles_subtnorm")
  profile_files <- list.files(profiles_subtnorm_dir,
                              pattern = "\\.(npz|parquet|csv)$")
  if (length(profile_files) == 0) {
    msg <- paste("❌ No profile files found in:", profile_main_dir)
    plc_log(msg, log_5, level = "❌ ERROR")
    stop(msg)
  }

  model_path <- file.path(system.file("model", package = "PRIME"), model_name)

  python_script_dir <- system.file("python", package = "PRIME")
  predict_script_path <- file.path(python_script_dir, "main.py")
  assertthat::assert_that(
    file.exists(predict_script_path),
    msg = paste("❌ Prediction script not found at:", predict_script_path)
  )

  assertthat::assert_that(
    file.exists(model_path),
    msg = paste("❌ Model file not found at:", model_path)
  )

  py_exec <- py_conf$python
  assertthat::assert_that(file.exists(py_exec),
                          msg = paste("❌ Python executable not found at:",
                                      py_exec))


  # Build Python command
  prediction_cmd <- c(
    py_exec, predict_script_path,
    "--script_dir", python_script_dir,
    "--profile_main_dir", profile_main_dir,
    "--combined_outdir", primeloci_tmp,
    "--model_path", model_path,
    "--log_file", log_5,
    "--name_prefix", name_prefix
  )

  if (!is.null(num_cores)) {
    prediction_cmd <- c(prediction_cmd, "--num_core", as.character(num_cores))
  }

  # Log the full Python command for debug before running
  plc_log(paste("🔧 Python command:",
                paste(shQuote(prediction_cmd), collapse = " ")),
          log_5, print_console = FALSE)

  plc_log("🔹 Running Python prediction script...", log_5)
  result <- tryCatch(
    {
      output <- system2(py_exec,
                        args = prediction_cmd[-1],
                        stdout = TRUE,
                        stderr = TRUE)
      attr(output, "status") <- 0
      output
    },
    error = function(e) {
      msg <- paste("❌ ERROR during prediction execution:", e$message)
      plc_log(msg, log_5, level = "❌ ERROR")
      attr(msg, "status") <- 1
      msg
    }
  )
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
    msg <- "❌ Prediction script failed. Check PRIMEloci-5.log for details."
    plc_log(msg, log_5, level = "❌ ERROR")
    stop(msg)
  }
  lapply(result, function(line) plc_log(paste("🔹", line), log_5, print_console = FALSE)) # nolint: line_length_linter.
  plc_log("✅ DONE :: Prediction step completed.", log_5)


  # _6_
  log_6 <- file.path(log_dir, "PRIMEloci-6.log")
  plc_log("\n\n\n 🚀 Running PRIMEloci -6 : Postprocessing prediction BEDs",
          log_6)

  bed_files <- find_bed_files_by_partial_name(primeloci_tmp,
                                              partial_name = postprocess_partial_name, # nolint: line_length_linter.
                                              log_file = log_6)
  if (length(bed_files) == 0) {
    msg <- paste("❌ No BED files found for postprocessing in", primeloci_tmp)
    plc_log(msg, log_6, level = "❌ ERROR")
    stop(msg)
  }

  plc_log(sprintf("📂 Found %d BED file(s) for processing.",
                  length(bed_files)),
          log_6)
  result_named_list <- lapply(seq_along(bed_files), function(i) {
    bed_file <- bed_files[i]
    basename_raw <- tools::file_path_sans_ext(basename(bed_file))
    pattern_match <- sub(paste0("^.*",
                                postprocess_partial_name,
                                "_(.*?)_combined.*$"),
                         "\\1",
                         basename_raw)

    sample_name <- if (identical(pattern_match, basename_raw)) {
      basename_raw
    } else {
      pattern_match
    }

    plc_log(paste("🔹 Processing:", bed_file), log_6)

    result_gr <- plc_trylog({
      coreovl_with_d(
        bed_file = bed_file,
        score_threshold = score_threshold,
        score_diff = score_diff,
        core_width = core_width,
        output_dir = outdir,
        num_cores = num_cores,
        log_file = log_6
      )
    }, log_file = log_6)

    if (!is.null(result_gr)) {
      return(list(name = sample_name, gr = result_gr))
    } else {
      plc_log(paste("⚠️ Skipped due to failure:", bed_file), log_6)
      return(NULL)
    }
  })


  # Filter out failed/null entries
  result_named_list <- Filter(Negate(is.null), result_named_list)

  if (length(result_named_list) == 0) {
    msg <- "❌ All postprocessing attempts failed or returned NULL."
    plc_log(msg, log_6, level = "❌ ERROR")
    stop(msg)
  }
  plc_log(sprintf("✅ DONE :: Postprocessed %d file(s) successfully.",
                  length(result_named_list)), log_6)

  # Final return object
  if (length(result_named_list) == 1) {
    result_gr_final <- result_named_list[[1]]$gr
  } else {
    result_gr_final <- GenomicRanges::GRangesList(
      setNames(
        lapply(result_named_list, `[[`, "gr"),
        vapply(result_named_list, `[[`, character(1), "name")
      )
    )
  }

  plc_log("✅ DONE :: Postprocessing completed.", log_6)
  message("✅✅✅ PRIMEloci pipeline completed successfully !!!!!")

  on.exit({
    if (!keep_tmp) {
      if (dir.exists(primeloci_tmp)) {
        unlink(primeloci_tmp, recursive = TRUE, force = TRUE)
        message(sprintf("Temporary directory '%s' has been cleaned up.",
                        primeloci_tmp))
      }
    }
  }, add = TRUE)

  return(result_gr_final)
}
