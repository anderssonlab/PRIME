#' Run the full plc pipeline for regulatory element prediction
#'
#' This function executes the complete PRIMEloci pipeline
#' to predict regulatory elements (enhancers and promoters)
#' from CTSS CAGE data. It integrates CAGE-derived tag clustering,
#' feature preparation, prediction using a pre-trained LightGBM model,
#' and post-processing to output high-confidence non-overlapping loci.
#'
#' The function is designed for users who
#' prefer a single-command workflow from a `RangedSummarizedExperiment` input
#' to the final genomic predictions.
#'
#' The pipeline was originally developed for human genome (hg38) CAGE data,
#' but can be adapted for other genomes with similar CAGE annotations.
#'
#' @param ctss_rse A `RangedSummarizedExperiment` object
#' containing CTSS-level expression data.
#' @param python_path Path to the Python environment to use.
#' This can be one of the following:
#' \itemize{
#'   \item A full path to a Python binary (e.g., `"/usr/bin/python3"`),
#'   \item A path to a virtualenv directory (must contain `bin/`),
#'   \item The name of a conda environment (e.g., `"prime-env"`),
#'   \item A full path to a conda environment directory.
#' }
#' The specified path must exist and be valid.
#' Default is `"~/.virtualenvs/prime-env"`.
#'
#' \strong{Important:} If Python is already initialized
#' (e.g., in RStudio or a long-running session),
#' changing the Python environment from within the function
#' will not take effect. To guarantee that the correct Python is used 
#' (especially when pointing to `"/usr/bin/python3"`),
#' set the environment variable `RETICULATE_PYTHON`
#' before starting R or RStudio.
#' Alternatively, call [configure_plc_python()] early in the session.
#'
#' @param score_threshold Minimum score threshold for core region predictions.
#' Must be between 0 and 1. Default is `0.75`.
#' @param score_diff Minimum score difference required between merged regions.
#' Must be non-negative and less than `score_threshold`. Default is `0.1`.
#' @param num_cores Number of cores to use for parallel processing.
#' Must be a positive integer or `NULL` to auto-detect. Default is `NULL`.
#' @param keep_tmp Logical. Whether to keep intermediate files
#' (e.g., profiles and temp folders). Default is `FALSE`.
#' @param log_dir Optional. Directory path
#' where a log file named `"PRIMEloci.log"`
#' will be written. If `NULL`, logs will be printed to the R console.
#' Default is `NULL`.
#'
#' @return A `GRanges` or `GRangesList` object containing
#' the final predicted loci after postprocessing.
#'
#' @details
#' The PRIMEloci pipeline includes the following steps:
#' \enumerate{
#'   \item \strong{Identifying Tag Clusters (TCs)}:
#' Identify tag clusters (TCs)
#' from the extracted CTSS data using the \pkg{CAGEfightR} package.
#'   \item \strong{Sliding Through TCs}:
#' Slide through the identified TCs (default window size = 20)
#' to create tiled regions for downstream processing.
#'   \item \strong{Creating Normalized Profiles}:
#' Generate normalized transcriptional profiles
#' suitable for input into the prediction model.
#'   \item \strong{Predicting Profile Probabilities}:
#' Use pre-trained PRIMEloci LightGBM models
#' to assign probabilities to each region,
#' indicating likelihood of being a regulatory element.
#'   \item \strong{Post-Processing}:
#' Refine and filter model predictions
#' using score thresholds and additional criteria.
#' Output non-overlapping core regulatory loci
#' in `GRanges` or `GRangesList` object, and (optional) BED file format
#' for further analysis.
#' }
#'
#' @section Python Environment:
#' This function attempts to configure Python using the `python_path` argument.
#' However, due to reticulate's behavior, Python must be configured
#' before initialization. If using a system Python path
#' (e.g., `"/usr/bin/python3"`), set `RETICULATE_PYTHON` before launching R.
#' For virtualenvs and conda environments,
#' configuration within the session usually works unless Python has already
#' been initialized.
#'
#' If `keep_tmp = FALSE`, temporary files will be removed
#' after the pipeline completes.
#'
#' @export
#'
#' @import GenomicRanges
#' @import assertthat
PRIMEloci <- function(
    ctss_rse,
    python_path = "~/.virtualenvs/prime-env",
    score_threshold = 0.75,
    score_diff = 0.1,
    num_cores = NULL,
    keep_tmp = FALSE,
    log_dir = NULL) {

  # Validate inputs

  # Check ctss_rse
  assertthat::assert_that(
    methods::is(ctss_rse, "RangedSummarizedExperiment"),
    msg = "`ctss_rse` must be a RangedSummarizedExperiment object."
  )

  # Set internal temporary output directory
  outdir <- file.path(tempdir(), "PRIMEloci_output")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message(sprintf("📁 Temporary output directory: %s", outdir))

  # Temporary working directories
  primeloci_tmp <- file.path(outdir, "PRIMEloci_tmp")
  dir.create(primeloci_tmp, recursive = TRUE, showWarnings = FALSE)

  # Set up logging
  if (is.null(log_dir)) {
    log_target <- stdout()  # Log to R console
  } else {
    log_dir <- normalizePath(path.expand(log_dir), mustWork = FALSE)

    assertthat::assert_that(
      is.character(log_dir),
      length(log_dir) == 1,
      msg = "`log_dir` must be a single character path or NULL."
    )

    if (!dir.exists(log_dir)) {
      dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    }

    log_target <- file.path(log_dir, "PRIMEloci.log")
    message(sprintf("📝 Log file will be saved to: %s", log_target))
  }

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

  py_conf <- configure_plc_python(python_path = python_path,
                                  log_target = log_target)

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

  plc_log("🚀 Starting PRIMEloci pipeline", log_target)

  # _2_
  plc_log("\n\n\n 🚀 Running PRIMEloci: get extended tc or validate the tc object provided", # nolint: line_length_linter.
          log_target)
  plc_log(sprintf("🕒 Pipeline started at: %s", Sys.time()), log_target)

  if (is.null(tc_object)) {
    plc_log("🔹 Creating tc object\n", log_target)
    tc_grl <- get_tcs_and_extend_fromthick(ctss_rse,
                                           ext_dis = ext_dis)
    if (save_tc) saveRDS(tc_grl, file = file.path(outdir, tc_object_name))
  } else {
    tc_grl <- tc_object
  }

  plc_log("🔹 Validating tc object\n", log_target)
  validate_tc <- validate_tc_object(tc_grl, ctss_rse, ext_dis = ext_dis)

  if (!validate_tc) {
    msg <- "\nTC object validation failed. Ensure the TC object is valid."
    plc_log(msg, log_target, level = "❌ ERROR")
    stop(msg)
  }
  plc_log("✅ DONE :: TC object is validated and ready to use.", log_target)


  # _3_
  plc_log("\n\n\n 🚀 Running PRIMEloci: sliding windows covering reduced TC regions", log_target) # nolint: line_length_linter.

  if (inherits(tc_grl, "GenomicRanges::GRanges")) {
    start_time <- Sys.time()
    tc_sliding_window_grl <- tc_sliding_window(tc_grl,
                                               sld_by = sld_by,
                                               ext_dis = ext_dis,
                                               log_file = log_target,
                                               num_cores = num_cores)
    plc_log(sprintf("⏱️ Time taken: %.2f minutes",
                    as.numeric(difftime(Sys.time(),
                                        start_time,
                                        units = "mins"))),
            log_target)
  } else if (inherits(tc_grl, "GenomicRanges::GRangesList") ||
               inherits(tc_grl, "CompressedGRangesList")) {
    tc_sliding_window_grl <- lapply(seq_along(tc_grl), function(i) {
      start_time <- Sys.time()
      gr_name <- if (!is.null(names(tc_grl))) names(tc_grl)[i] else paste0("Sample_", i) # nolint: line_length_linter.
      plc_log(sprintf("🔹 Processing: %s", gr_name), log_target)
      result <- tc_sliding_window(tc_grl[[i]],
                                  sld_by = sld_by,
                                  ext_dis = ext_dis,
                                  log_file = log_target,
                                  num_cores = num_cores)
      plc_log(sprintf("⏱️ Time taken: %.2f minutes",
                      as.numeric(difftime(Sys.time(),
                                          start_time,
                                          units = "mins"))),
              log_target)
      result
    })

    if (!is.null(tc_sliding_window_grl) &&
          length(tc_sliding_window_grl) > 0) {
      tc_sliding_window_grl <- GenomicRanges::GRangesList(tc_sliding_window_grl) # nolint: line_length_linter.
      names(tc_sliding_window_grl) <- names(tc_grl)
    } else {
      msg <- "Processed TC object list is empty. Ensure tc_grl contains valid data." # nolint: line_length_linter.
      plc_log(msg, log_target, level = "❌ ERROR")
      stop(msg)
    }
  } else {
    msg <- "tc_grl must be either a GRanges, GRangesList, or CompressedGRangesList object." # nolint: line_length_linter.
    plc_log(msg, log_target, level = "❌ ERROR")
    stop(msg)
  }

  if (save_sld_tc) {
    plc_log(sprintf("Saving TC objects to PRIMEloci_tmp .."),
            log_target)
    saveRDS(tc_sliding_window_grl,
            file.path(primeloci_tmp, sld_object_name))
  }
  plc_log("✅ DONE :: Sliding window TC object is saved to PRIMEloci_tmp.",
          log_target)
  tc_for_profile <- tc_sliding_window_grl

  assertthat::assert_that(!is.null(tc_for_profile),
                          msg = "❌ tc_for_profile is NULL at return. Something failed.") # nolint: line_length_linter.


  # _4_
  plc_log("\n\n\n 🚀 Running PRIMEloci: compute count & normalized profiles for each sample", # nolint: line_length_linter.
          log_target)

  plc_profile(
    ctss_rse,
    tc_for_profile,
    outdir,
    profile_dir_name,
    file_type = file_type,
    python_path = py_conf$python,
    addtn_to_filename = addtn_to_filename,
    save_count_profiles = save_count_profiles,
    num_cores = num_cores,
    log_file = log_target,
    ext_dis
  )


  # _5_
  plc_log("\n\n\n 🚀 Running PRIMEloci: Prediction using PRIMEloci model",
          log_target)

  profile_main_dir <- file.path(primeloci_tmp,
                                profile_dir_name)

  profiles_subtnorm_dir <- file.path(profile_main_dir, "profiles_subtnorm")
  profile_files <- list.files(profiles_subtnorm_dir,
                              pattern = "\\.(npz|parquet|csv)$")
  if (length(profile_files) == 0) {
    msg <- paste("❌ No profile files found in:", profile_main_dir)
    plc_log(msg, log_target, level = "❌ ERROR")
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
    "--log_file", log_target,
    "--name_prefix", name_prefix
  )

  if (!is.null(num_cores)) {
    prediction_cmd <- c(prediction_cmd, "--num_core", as.character(num_cores))
  }

  # Log the full Python command for debug before running
  plc_log(paste("🔧 Python command:",
                paste(shQuote(prediction_cmd), collapse = " ")),
          log_target, print_console = FALSE)

  plc_log("🔹 Running Python prediction script...", log_target)
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
      plc_log(msg, log_target, level = "❌ ERROR")
      attr(msg, "status") <- 1
      msg
    }
  )
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
    msg <- "❌ Prediction script failed. Check PRIMEloci.log for details."
    plc_log(msg, log_target, level = "❌ ERROR")
    stop(msg)
  } else {
    plc_log("✅ DONE :: Prediction script executed successfully.", log_target)
  }


  # _6_
  plc_log("\n\n\n 🚀 Running PRIMEloci: Postprocessing prediction BEDs",
          log_target)

  bed_files <- find_bed_files_by_partial_name(primeloci_tmp,
                                              partial_name = postprocess_partial_name, # nolint: line_length_linter.
                                              log_file = log_target)
  if (length(bed_files) == 0) {
    msg <- paste("❌ No BED files found for postprocessing in", primeloci_tmp)
    plc_log(msg, log_target, level = "❌ ERROR")
    stop(msg)
  }

  plc_log(sprintf("📂 Found %d BED file(s) for processing.",
                  length(bed_files)),
          log_target)
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

    result_gr <- coreovl_with_d(bed_file = bed_file,
                                score_threshold = score_threshold,
                                score_diff = score_diff,
                                core_width = core_width,
                                return_gr = TRUE,
                                output_dir = NULL,
                                num_cores = num_cores,
                                log_file = log_target)

    if (!is.null(result_gr)) {
      list(name = sample_name, gr = result_gr)
    } else {
      plc_log(paste("⚠️ Skipped due to failure:", bed_file), log_target)
      NULL
    }
  })

  # Filter out failed/null entries
  result_named_list <- Filter(Negate(is.null), result_named_list)

  if (length(result_named_list) == 0) {
    msg <- "❌ All postprocessing attempts failed or returned NULL."
    plc_log(msg, log_target, level = "❌ ERROR")
    stop(msg)
  }
  plc_log(sprintf("✅ DONE :: Postprocessed %d file(s) successfully.",
                  length(result_named_list)), log_target)

  sample_names <- disambiguate_sample_names(result_named_list, log_target)

  # Final return object
  if (length(result_named_list) == 1) {
    result_gr_final <- result_named_list[[1]]$gr
  } else {
    result_gr_final <- GenomicRanges::GRangesList(
      setNames(
        lapply(result_named_list, `[[`, "gr"),
        sample_names
      )
    )
  }
  
  on.exit({
    if (!keep_tmp) {
      if (dir.exists(outdir)) {
        unlink(outdir, recursive = TRUE, force = TRUE)
        message(sprintf("Temporary directory '%s' has been cleaned up.",
                        outdir))
      }
    }
  }, add = TRUE)

  plc_log("✅ DONE :: Postprocessing completed.", log_target)
  message("✅✅✅ PRIMEloci pipeline completed successfully !!!!!")
  plc_log(sprintf("🏁 Pipeline completed at: %s", Sys.time()), log_target)

  return(result_gr_final)
}
