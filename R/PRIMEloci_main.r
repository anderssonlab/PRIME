#' Run the PRIMEloci Workflow
#'
#' This function encapsulates the entire PRIMEloci workflow,
#' including creating sliding windows, generating profiles,
#' running a Python script for profile prediction,
#' and processing the results into a `GRangesList`.
#'
#' @param ctss_rse A `SummarizedExperiment` object containing CTSS data.
#' @param config_file Character. Path to the configuration file in YAML format.
#' Default is "config_R_PRIMEloci.yaml".
#'
#' @return A `GRangesList` containing the processed results.
#'
#' @importFrom yaml read_yaml
#' @importFrom GenomicRanges GRangesList sort
#' @importFrom future.apply future_lapply
#' @importFrom parallel detectCores
#' @importFrom reticulate py_run_string
#' @importFrom argparse ArgumentParser
#' @export
PRIMEloci <- function(ctss_rse, config_file = "config_R_PRIMEloci.yaml") {

  ### **STEP 1: Set Up Directories & Load Configuration**
  tmp_dir <- tempdir()
  dir.create(tmp_dir)

  config <- list(
    output_dir = tmp_dir,
    profile_main_dir = "profiles",
    python_script_dir = system.file("python", package = "PRIME"),
    model_path = system.file("models", "PRIMEloci_GM12878_wt10M.sav", package = "PRIME"),
    prefix_out_name = "modelPred",
    profile_sub_dir = "tcs",
    profile_file_type = "parquet",
    threshold = 0.2,
    save_count_profiles = TRUE,
    ext_dis = 200,
    partial_name = "pred_slt.*\\.bed"
  )

  # Prepare output directories
  prediction_dir <- file.path(config$output_dir, config$profile_main_dir, "predictions", config$profile_sub_dir)
  outdir_main_name <- c("metadata", "profiles", "profiles_subtnorm", "predictions")
  prep_profile_dir(config$output_dir, config$profile_main_dir, outdir_main_name, config$profile_sub_dir)

  ### **STEP 2: Extract Tag Clusters & Extend**
  writeLines("\nExtracting and extending TCs..")
  tc_grl <- suppressWarnings(get_tcs_and_extend_fromthick(ctss_rse, ext_dis=config$ext_dis))

  ### **STEP 3: Sliding Window Processing**
  writeLines("\nApplying Sliding Window to Tag Clusters..")
  tc_sliding_window_grl <- lapply(seq_along(tc_grl), function(i) {
    print(paste("Processing:", names(tc_grl)[i]))
    tc_sliding_window(tc_grl[[i]], slide_by = 20, expand_by = 200, num_cores = 2)
  })
  tc_sliding_window_grl <- GenomicRanges::GRangesList(tc_sliding_window_grl)

  ### **STEP 4: Generate Profiles**
  writeLines("\nGenerating Profiles..")
  PRIMEloci_profile_2(
    ctss_rse,
    tc_sliding_window_grl,
    config$output_dir,
    config$profile_main_dir,
    config$profile_sub_dir,
    config$ext_dis,
    save_count_profiles = config$save_count_profiles,
    file_type = config$profile_file_type
  )

  ### **STEP 5: Run Python Prediction Script**
  writeLines("\nRunning Python Profile Prediction..")
  run_prediction_python_script(
    script_path = system.file("python", "_predict_profile_probabilities.py", package = "PRIME"),
    script_dir = config$python_script_dir,
    profile_main_dir = file.path(config$output_dir, config$profile_main_dir),
    profile_sub_dir = config$profile_sub_dir,
    model_path = config$model_path,
    name_prefix = config$prefix_out_name,
    threshold = config$threshold,
    file_format = config$profile_file_type
  )

  ### **STEP 6: Postprocessing (Core Overlap & Refinement)**
  writeLines("\nPostprocessing Predicted Profiles..")
  bed_file <- list.files(prediction_dir, pattern = config$partial_name, full.names = TRUE)
  processed_gr <- process_all_files(bed_file)

  # Remove temporary files
  unlink(tmp_dir, recursive = TRUE)

  writeLines("PRIMEloci Workflow Completed!")

  return(processed_gr)
}
