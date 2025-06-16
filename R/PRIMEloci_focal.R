#' Run the PRIMEloci focal Pipeline on CTSS and Identified Regions
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
#'   in the desired environment. Default is NULL.
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
PRIMEloci_focal <- function(
    ctss_rse,
    tc_gr,
    python_path = NULL,
    num_cores = NULL,
    keep_tmp = FALSE,
    log_dir = NULL) {


  # setting
  profile_dir_name <- "PRIMEloci_profiles"
  postprocess_partial_name <- "pred_all"
  save_count_profiles <- FALSE
  ext_dis <- 200
  file_type <- "npz"
  addtn_to_filename <- ""
  name_prefix <- "PRIMEloci"
  model_name <- "PRIMEloci_GM12878_model_1.0.sav"


  # Validate inputs

  # Check ctss_rse
  assertthat::assert_that(
    methods::is(ctss_rse, "RangedSummarizedExperiment"),
    msg = "❌ `ctss_rse` must be a RangedSummarizedExperiment object."
  )

  # Set internal temporary output directory
  outdir <- file.path(tempdir())
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  plc_message(sprintf("📁 Temporary output directory: %s", outdir))

  # Logging setup
  if (is.null(log_dir)) {
    log_target <- stdout()  # Console-only
  } else {
    log_dir <- normalizePath(path.expand(log_dir), mustWork = FALSE)

    assertthat::assert_that(
      is.character(log_dir),
      length(log_dir) == 1,
      msg = "❌ `log_dir` must be a single character path or NULL."
    )

    if (!dir.exists(log_dir)) {
      dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    }

    log_target <- file.path(log_dir, "PRIMEloci.log")
    plc_message(sprintf("📝 Log file will be saved to: %s", log_target))

    # This captures stdout (e.g., cat, print, plc_log) to both console + file
    sink(log_target, append = TRUE, split = TRUE)

    on.exit({
      sink(NULL)  # close stdout sink
    }, add = TRUE)
  }

  if (!is.null(num_cores)) {
    assertthat::assert_that(
      is.numeric(num_cores),
      num_cores %% 1 == 0,
      num_cores > 0,
      msg = "❌ `num_cores` must be a positive integer or NULL."
    )
  }

  if (is.null(num_cores)) {
    num_cores <- max(1, min(25, parallel::detectCores() %/% 2))
  }
  if (num_cores == 1) {
    processing_method <- "callr"
    plc_message("⚠️ num_workers was set to 1. Using callr backend: tasks will run sequentially (despite using multiple R sessions).") # nolint: line_length_linter.
  } else {
    processing_method <- plc_detect_parallel_plan()
  }

  plc_message("\n")
  plc_message("🚀 Setting up Python environment")
  if (is.null(python_path)) {
    py <- reticulate::import("sys")
    python_path <- py$executable
  }
  py_conf <- plc_configure_python(python_path = python_path)
  check_npz <- plc_test_scipy_save_npz()
  if (!check_npz) {
    plc_message("⚠️ Falling back to .parquet format")
    file_type <- "parquet"
  } else {
    plc_message("✅ Using .npz format")
    file_type <- "npz"
  }

  plc_message("\n")
  plc_message("🚀 Starting PRIMEloci focal pipeline")

  # _2_
  plc_message("\n")
  plc_message("🚀 Validating the tc granges object provided")

  len_vec <- ext_dis * 2 + 1

  # Check if tc_object is a GRanges
  assertthat::assert_that(
    inherits(tc_gr, "GRanges"),
    msg = "❌ The object must be a GRanges object."
  )

  # Ensure all regions have the correct width
  if (all(GenomicRanges::width(tc_gr) != len_vec)) {
    msg <- paste("⚠️ All regions in the object (GRanges) must have width",
                 len_vec,
                 " : extend 401 bp from thick if existed")
    tc_gr <- plc_extend_fromthick(tc_gr = tc_gr,
                                  ext_dis = ext_dis)
  } else {
    msg <- paste("✅ All regions in the object (GRanges) have width", len_vec) # nolint: line_length_linter.
  }
  plc_message(msg)

  validate_tc <- plc_validate_tc_object(tc_gr, ctss_rse, ext_dis = ext_dis)

  if (!validate_tc) {
    plc_error("❌ TC object validation failed. Ensure the TC object is valid.")
  }
  plc_message("✅ DONE :: TC object is validated and ready to use.")


  # _4_
  plc_message("\n")
  plc_message("🚀 Computing count & normalized profiles for each sample")

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
    processing_method = processing_method,
    ext_dis
  )


  # _5_
  plc_message("\n")
  plc_message("🚀 Predicting probability using PRIMEloci model")

  profile_main_dir <- file.path(outdir,
                                profile_dir_name)
  profiles_subtnorm_dir <- file.path(profile_main_dir, "profiles_subtnorm")
  profile_files <- list.files(profiles_subtnorm_dir,
                              pattern = "\\.(npz|parquet|csv)$")
  if (length(profile_files) == 0) {
    plc_error(paste("❌ No profile files found in:", profile_main_dir))
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
    "--combined_outdir", outdir,
    "--model_path", model_path,
    "--log_file", log_target,
    "--name_prefix", name_prefix
  )

  if (!is.null(num_cores)) {
    prediction_cmd <- c(prediction_cmd, "--num_core", as.character(num_cores))
  }

  # Log the full Python command for debug before running
  plc_message(paste("🔧 Python command:",
                    paste(shQuote(prediction_cmd), collapse = " ")))

  plc_message("🔹 Running Python prediction script...")
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
      msg <- paste("❌ ERROR during prediction execution: ", e$message)
      plc_message(msg)
      attr(msg, "status") <- 1
      msg
    }
  )
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
    plc_error("❌ Prediction script failed. Check PRIMEloci.log for details.")
  } else {
    plc_message("✅ DONE :: Prediction script executed successfully.")
  }

  plc_message("\n")
  plc_message("🚀 Importing prediction BEDs")

  final_rse <- plc_focal_prediction_to_rse(outdir,
                                           postprocess_partial_name = postprocess_partial_name) # nolint: line_length_linter.

  on.exit({
    if (!keep_tmp) {
      if (dir.exists(outdir)) {
        unlink(outdir, recursive = TRUE, force = TRUE)
        plc_message((sprintf("🧹 Temporary directory '%s' has been cleaned up.",
                             outdir)))
      }
    }
  }, add = TRUE)

  plc_message("\n")
  plc_message("✅✅✅ PRIMEloci focal pipeline completed successfully !!!!!")
  plc_message(sprintf("🏁 Pipeline completed at: %s", Sys.time()))
  plc_message("\n")

  return(final_rse)
}
