#' Run PRIMEloci on ctss_rse_chr16to17.rds from extdata
#'
#' @param outdir Output directory for the pipeline
#' @param python_path Path to Python virtual environment
#' @param num_cores Number of CPU cores to use
#'
#' @return GRanges or GRangesList
#' @export
run_PRIMEloci_example <- function(python_path = NULL,
                                  log_dir = NULL,
                                  keep_tmp = FALSE) {
  rds_path <- system.file("extdata",
                          "ctss_rse_chr16to17.rds",
                          package = "PRIME")
  stopifnot(file.exists(rds_path))

  ctss_rse <- readRDS(rds_path)

  result <- PRIMEloci(ctss_rse = ctss_rse,
                      python_path = python_path,
                      log_dir = log_dir,
                      keep_tmp = keep_tmp)

  return(result)
}
