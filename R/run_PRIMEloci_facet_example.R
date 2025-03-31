#' Run an Example of the PRIMEloci Facet Pipeline
#'
#' This function demonstrates how to run the `PRIMEloci_facet()` pipeline
#' using example CTSS and region data bundled with the `PRIME` package.
#' It loads pre-packaged test data, runs the pipeline, and returns the result.
#'
#' @param python_path Character path to the Python binary
#'   in the desired environment. Default is `"~/.virtualenvs/prime-env"`.
#' @param log_dir Optional character path to save the log file.
#'   Default is `"/Users/natsudanav/Desktop/PRIMEloci_facet.log"`.
#' @param keep_tmp Logical. If `TRUE`, temporary files and folders are retained.
#'   Default is `FALSE`.
#'
#' @return A `GRanges` object or a `GRangesList` object,
#'   depending on the number of samples in the example dataset.

#' @export
run_PRIMEloci_facet_example <- function(python_path = "~/.virtualenvs/prime-env",
                                        log_dir = NULL,
                                        keep_tmp = FALSE) {

  rds_ctss <- system.file("extdata",
                          "ctss_rse_chr16to17.rds",
                          package = "PRIME")
  rds_tc <- system.file("extdata",
                        "predicted_regions_gr.rds",
                        package = "PRIME")

  stopifnot(file.exists(rds_ctss))
  stopifnot(file.exists(rds_tc))

  ctss_rse <- readRDS(rds_ctss)
  tc_gr <- readRDS(rds_tc)

  result <- PRIMEloci_facet(ctss_rse = ctss_rse,
                            tc_gr = tc_gr,
                            python_path = python_path,
                            log_dir = log_dir,
                            keep_tmp = keep_tmp)

  return(result)
}