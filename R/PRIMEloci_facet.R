#' Run the PRIMEloci Facet Pipeline on CTSS and Identified Regions
#'
#' This function executes a partial PRIMEloci pipeline,
#' which includes the steps of validating and extending existing regions—
#' such as those predicted from pooled PRIMEloci data—
#' followed by profile generation, model prediction, and BED result import.
#'
#' @param ctss_rse A `RangedSummarizedExperiment` object representing CTSS data.
#' @param tc_gr A `GRanges` object representing identified regions.
#'   Regions will be extended to 401 bp width if needed.
#' @param python_path Character path to the Python binary
#'   in the desired environment. Default is `"~/.virtualenvs/prime-env"`.
#' @param num_cores Optional integer specifying the number of CPU cores
#'   to use for parallel steps.
#' @param keep_tmp Logical. If `TRUE`, temporary files and folders
#' will be retained.
#'   Default is `FALSE`.
#' @param log_dir Optional path to save a log file.
#'   If `NULL`, logs are printed to the console.
#'
#' @return A `GRanges` object if one sample was processed,
#'   or a `GRangesList` object if multiple samples were processed.
#'
#' @details
#' This function supports downstream analysis by:
#' - Validating and, if necessary, extending input regions to 401 bp
#' - Generating normalized CTSS-based profiles per sample
#' - Running predictions using a pre-trained PRIMEloci model via Python
#' - Importing the prediction results from BED files
#'
#' Temporary files are stored in a subdirectory of `tempdir()`
#' and removed unless `keep_tmp = TRUE`.
#'
#' @import GenomicRanges
#' @import assertthat
#' @export
PRIMEloci_facet <- function(
    ctss_rse,
    tc_gr,
    python_path = "~/.virtualenvs/prime-env",
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
  profile_dir_name <- "PRIMEloci_profiles"
  postprocess_partial_name <- "pred_all"
  save_count_profiles <- FALSE
  ext_dis <- 200
  file_type <- "npz"
  addtn_to_filename <- ""
  name_prefix <- "PRIMEloci"
  model_name <- "PRIMEloci_GM12878_model_1.0.sav"

  plc_log("🚀 Starting PRIMEloci facet pipeline", log_target)

  # _2_
  plc_log("\n\n\n 🚀 Validating the tc granges object provided",
          log_target)

  len_vec <- ext_dis * 2 + 1

  # Check if tc_object is a GRanges
  assertthat::assert_that(
    inherits(tc_gr, "GRanges"),
    msg = "\n❌ the object must be a GRanges object"
  )

  # Ensure all regions have the correct width
  if (all(GenomicRanges::width(tc_gr) != len_vec)) {
    msg <- paste("\n⚠️ All regions in the object (GRanges) must have width",
                 len_vec,
                 " : extend 401 bp from thick if existed")
    tc_gr <- extend_fromthick(tc_gr = tc_gr,
                              ext_dis = ext_dis)
  } else {
    msg <- paste("\n✅ All regions in the object (GRanges) have width", len_vec) # nolint: line_length_linter.
  }
  plc_log(msg, log_target)

  validate_tc <- validate_tc_object(tc_gr, ctss_rse, ext_dis = ext_dis)

  if (!validate_tc) {
    msg <- "\nTC object validation failed. Ensure the TC object is valid."
    plc_log(msg, log_target, level = "❌ ERROR")
    stop(msg)
  }
  plc_log("✅ DONE :: TC object is validated and ready to use.", log_target)


  # _4_
  plc_log("\n\n\n 🚀 Computing count & normalized profiles for each sample", # nolint: line_length_linter.
          log_target)

  plc_profile(
    ctss_rse,
    tc_gr,
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
  plc_log("\n\n\n 🚀 Predicting probability using PRIMEloci model",
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

  plc_log("\n\n\n 🚀 Importing prediction BEDs",
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
    # Load and prepare data
    bed <- load_bed_file(bed_file)
    gr <- create_granges_from_bed(bed)

    sample_name <- if (identical(pattern_match, basename_raw)) {
      basename_raw
    } else {
      pattern_match
    }

    if (!is.null(gr)) {
      list(name = sample_name, gr = gr)
    } else {
      plc_log(paste("⚠️ Skipped due to failure:", bed_file), log_target)
      NULL
    }
  })

  # Filter out failed/null entries
  result_named_list <- Filter(Negate(is.null), result_named_list)

  if (length(result_named_list) == 0) {
    msg <- "❌ All attempts failed or returned NULL."
    plc_log(msg, log_target, level = "❌ ERROR")
    stop(msg)
  }
  plc_log(sprintf("✅ DONE :: importing %d file(s) successfully.",
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
        plc_log(sprintf("🧹 Temporary directory '%s' has been cleaned up.",
                        outdir),
                log_target)
      }
    }
  }, add = TRUE)
  
  message("✅✅✅ PRIMEloci facet pipeline completed successfully !!!!!")
  plc_log(sprintf("🏁 Pipeline completed at: %s", Sys.time()), log_target)

  return(result_gr_final)
}
