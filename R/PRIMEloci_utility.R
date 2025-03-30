#' Log a timestamped message with level (and optional console output)
#'
#' Writes a log entry to a file with a timestamp and log level.
#' Optionally prints the same message to the R console.
#'
#' @param message A character string to log.
#' @param log_file Path to the log file.
#' @param level Log level string (e.g., "INFO", "WARN", "ERROR").
#' Default is "INFO".
#' @param print_console Logical. If TRUE, also prints to console
#' (default: TRUE).
#'
#' @return None. Writes to the specified log file.
plc_log <- function(message,
                    log_file,
                    level = "INFO",
                    print_console = FALSE) {
  timestamp <- base::format(base::Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- base::sprintf("[%s] %s %s", level, timestamp, message)
  base::cat(line, "\n", file = log_file, append = TRUE)
  if (isTRUE(print_console)) base::message(line)
}

#' Capture and log messages, warnings, and errors from a code block
#'
#' @param expr Code to evaluate.
#' @param log_file Path to the log file.
#' @param print_console Logical. Whether to also print output to console.
#'
#' @return The result of the evaluated expression, or NULL on error.
plc_trylog <- function(expr, log_file, print_console = TRUE) {
  withCallingHandlers(
    tryCatch(
      {
        eval(expr)
      },
      error = function(e) {
        plc_log(conditionMessage(e),
                log_file,
                level = "❌ ERROR",
                print_console = print_console)
        NULL  # return NULL on error
      }
    ),
    warning = function(w) {
      plc_log(conditionMessage(w),
              log_file,
              level = "⚠️ WARN",
              print_console = print_console)
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      plc_log(conditionMessage(m),
              log_file,
              level = "INFO",
              print_console = print_console)
      invokeRestart("muffleMessage")
    }
  )
}

#' Internal utility to check that all required Python packages are installed
#' and meet the minimum version requirements. Raises an error if any package
#' is missing or incompatible. Not exported.
check_python_dependencies <- function(required_packages) {
  message("🔍 Checking Python dependencies...\n")
  for (pkg in names(required_packages)) {
    min_version <- required_packages[[pkg]]
    result <- tryCatch({
      mod <- reticulate::import(pkg, delay_load = FALSE)
      version <- as.character(mod$`__version__`)
      version_msg <- paste0("found version ", version)
      if (!is.null(min_version) &&
            utils::compareVersion(version, min_version) < 0) {
        stop(paste0("⚠️  Package '", pkg, "' ", version_msg,
                    " < required ", min_version))
      }
      message(paste0("✅ ", pkg, ": ", version_msg))
      TRUE
    }, error = function(e) {
      message(paste0("❌ ",
                     pkg,
                     ": not found or import failed (", e$message, ")"))
      FALSE
    })
    if (!result) {
      stop(paste0("🛑 Missing or incompatible Python package: '", pkg, "'"))
    }
  }
  message("\n✅ All required Python packages are available.\n")
}

#' Internal test to confirm that scipy.sparse.save_npz() works via reticulate.
#' This ensures compatibility between R, reticulate, scipy, and numpy.
plc_test_scipy_save_npz <- function(log_target = NULL) {
  tryCatch({
    test_matrix <- Matrix::rsparsematrix(10, 10, density = 0.1)
    scipy <- reticulate::import("scipy.sparse", delay_load = FALSE)

    if (inherits(test_matrix, "dgCMatrix")) {
      test_py <- scipy$csc_matrix(reticulate::tuple(test_matrix@x,
                                                    test_matrix@i,
                                                    test_matrix@p),
                                  shape = reticulate::tuple(test_matrix@Dim[1],
                                                            test_matrix@Dim[2]))
    } else if (inherits(test_matrix, "dgRMatrix")) {
      test_py <- scipy$csr_matrix(reticulate::tuple(test_matrix@x,
                                                    test_matrix@j,
                                                    test_matrix@p),
                                  shape = reticulate::tuple(test_matrix@Dim[1],
                                                            test_matrix@Dim[2]))
    } else {
      stop("Unsupported test matrix type for SciPy save check.")
    }

    tmpfile <- tempfile(fileext = ".npz")
    scipy$save_npz(tmpfile, test_py)
    unlink(tmpfile)
    plc_log("✅ SciPy sparse matrix save test passed.", log_target)
  }, error = function(e) {
    msg <- paste("❌ Critical error: scipy.sparse.save_npz() failed — cannot proceed.\n", # nolint: line_length_linter.
                 "➡️ Set RETICULATE_PYTHON before starting R if using system Python.\n", # nolint: line_length_linter.
                 "Full error:\n", e$message)
    plc_log(msg, log_target, level = "❌ ERROR")
    stop(msg)
  })
}

#' Configure Python environment for PRIMEloci
#'
#' Sets and verifies the Python environment
#' to be used in the PRIMEloci pipeline.
#' This includes specifying a path to a Python binary,
#' virtual environment, or conda environment,
#' checking for required Python packages and their versions,
#' and ensuring runtime compatibility.
#'
#' @details
#' In addition to package version checks,
#' this function performs a runtime test to verify
#' that `scipy.sparse.save_npz()` works correctly through `reticulate`.
#' This prevents downstream failures during profile saving
#' if Python is misconfigured, especially when
#' using system Python without setting `RETICULATE_PYTHON` early in the session.
#'
#' @param python_path Path to Python binary, virtualenv, or conda environment.
#' @param log_target Path to a log file for output (optional).
#'
#' @return Python configuration object from `reticulate::py_config()`
#' @export
configure_plc_python <- function(python_path = "~/.virtualenvs/prime-env",
                                       log_target = NULL) {
  python_path <- path.expand(python_path)

  if (reticulate::py_available(initialize = FALSE)) {
    warning("Python has already been initialized — changes may not take effect. Set RETICULATE_PYTHON before starting R or call this early in the session.") # nolint: line_length_linter.
  }

  if (!is.null(python_path) && file.exists(python_path)) {
    if (grepl("bin/python", python_path) ||
          grepl("python[0-9.]*$", python_path)) {
      reticulate::use_python(python_path, required = TRUE)
    } else if (dir.exists(file.path(python_path, "bin"))) {
      reticulate::use_virtualenv(python_path, required = TRUE)
    } else if (basename(python_path) %in% reticulate::conda_list()$name) {
      reticulate::use_condaenv(python_path, required = TRUE)
    } else if (dir.exists(python_path) && grepl("conda", python_path)) {
      reticulate::use_condaenv(python_path, required = TRUE)
    } else {
      stop("Invalid python_path: Not recognized as a Python binary, virtualenv, or conda environment.") # nolint: line_length_linter.
    }
  } else {
    stop("python_path does not exist: Please provide a valid path to a Python binary, virtualenv, or conda environment.") # nolint: line_length_linter.
  }

  py_conf <- reticulate::py_config()

  # Logging
  if (!is.null(log_target)) {
    plc_log(sprintf("🔧 R version: %s", R.version.string), log_target)
    plc_log(sprintf("📦 Loaded Python: %s", py_conf$python), log_target)
  }

  # Required packages (match your validated env)
  required_packages <- list(
    numpy         = "2.0.2",
    scipy         = "1.13.1",
    pandas        = "2.2.3",
    joblib        = "1.4.2",
    sklearn       = "1.4.2",
    lightgbm      = "4.6.0",
    pyarrow       = "17.0.0",
    fastparquet   = "2024.11.0"
  )
  check_python_dependencies(required_packages)
  plc_test_scipy_save_npz(log_target)

  return(py_conf)
}

#' Get tag clusters and extend from thick positions
#'
#' This function calculates tag clusters from
#' a RangedSummarizedExperiment object \code{ctss_rse}
#' and extends them by a specified distance \code{ext_dis}
#' around the thick positions.
#'
#' @param ctss_rse A RangedSummarizedExperiment object containing CAGE data.
#' @param ext_dis An integer specifying the distance
#' to extend around thick positions (default is 200).
#'
#' @return A GenomicRanges::GRangesList object where each element corresponds to
#'         a tag cluster extended by \code{ext_dis} around the thick positions.
#'
#' @details
#' The function iterates over columns of \code{ctss_rse},
#' calculates pooled counts,
#' subsets clusters based on a score threshold (> 0),
#' and clusters them unidirectionally.
#' It then extends each cluster by \code{ext_dis} base pairs
#' around the thick positions
#' and stores the results in a GRangesList object.
#'
#' @examples
#' # Example usage with a RangedSummarizedExperiment object
#' # ctss_rse <- ...  # Load or create your RangedSummarizedExperiment object
#' # ext_dis <- 200   # Define your extension distance
#' # result <- get_tagclusters_and_extend_fromthick(ctss_rse, ext_dis)
#'
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @import IRanges
#' @import CAGEfightR
#' @import assertthat
#'
#' @export
get_tcs_and_extend_fromthick <- function(ctss_rse, ext_dis = 200) {

  # Assert that ctss_rse is a RangedSummarizedExperiment object
  assertthat::assert_that(
    inherits(ctss_rse, "RangedSummarizedExperiment"),
    msg = "ctss_rse must be a RangedSummarizedExperiment object."
  )

  # Get column names
  col_ctss_rse <- colnames(ctss_rse)

  # Initialize GRangesList
  tc_grl <- GenomicRanges::GRangesList()

  # Loop over column names
  for (i in col_ctss_rse) {

    writeLines(paste("Processing: ", i))

    # Extract data for current column
    ctss <- ctss_rse[, i]

    # Calculate pooled counts
    ctss <- CAGEfightR::calcPooled(ctss, inputAssay = "counts")

    # Subset by score > 0 (score column is created by calcPooled())
    ctss <- base::subset(ctss, score > 0)

    # Cluster unidirectionally
    object <- CAGEfightR::clusterUnidirectionally(ctss)

    # Create new ranges around thick positions
    new_ranges <- IRanges::IRanges(start = start(object$thick) - ext_dis,
                                   end = end(object$thick) + ext_dis)

    # Assign new ranges to object
    new_object <- object
    IRanges::ranges(new_object) <- new_ranges

    # Trim new object
    writeLines("Trimming out-of-bound ranges...")
    new_object <- GenomicRanges::trim(new_object)

    writeLines("Keep only prefered width...\n")
    len_vec <- ext_dis * 2 + 1

    new_object_widths <- GenomicRanges::width(new_object)
    new_object <- new_object[new_object_widths == len_vec]

    # Store in GRangesList
    tc_grl[[i]] <- new_object
  }

  return(tc_grl)
}

#' Validate a TC object
#'
#' This function ensures that a given `GenomicRanges::GRanges` or
#' `GenomicRanges::GRangesList` object meets the expected criteria.
#'
#' @param tc_object A `GenomicRanges::GRanges` or
#' `GenomicRanges::GRangesList` object.
#' @param ctss_rse The original CTSS dataset (for length and naming validation).
#' @param ext_dis An integer value for extension distance. Default is 200.
#'
#' @return TRUE if validation passes; otherwise, throws an error.
#'
#' #' @examples
#' # Validate an existing TC object:
#' validate_tc_object(tc_object, ctss_rse)
#'
#' @import GenomicRanges
#' @import assertthat
#'
#' @export
validate_tc_object <- function(tc_object, ctss_rse, ext_dis = 200) {

  len_vec <- ext_dis * 2 + 1

  # Check if tc_object is a valid GRanges or GRangesList
  assertthat::assert_that(
    inherits(tc_object, "GRanges") ||
      inherits(tc_object, "GRangesList") ||
      inherits(tc_object, "CompressedGRangesList"),
    msg = "\n❌ tc_object must be a GRanges, GRangesList, or CompressedGRangesList object" # nolint: line_length_linter.
  )

  if (inherits(tc_object, "GRanges")) {

    # Ensure all regions have the correct width
    assertthat::assert_that(
      all(GenomicRanges::width(tc_object) == len_vec),
      msg = paste("\n❌ All regions in tc_object (GRanges) must have width", len_vec) # nolint: line_length_linter.
    )

    # Check seqlevels without modifying
    sl_tc <- GenomeInfoDb::seqlevels(tc_object)
    sl_ctss <- GenomeInfoDb::seqlevels(ctss_rse)
    if (!setequal(sl_tc, sl_ctss)) {
      warning("⚠️ seqlevels differ between tc_object (GRanges) and ctss_rse — pruning is NOT applied during validation.") # nolint: line_length_linter.
      message("  → tc_object seqlevels: ", paste(sl_tc, collapse = ", "))
      message("  → ctss_rse seqlevels:  ", paste(sl_ctss, collapse = ", "))
    } else {
      message("seqlevels match between tc_object and ctss_rse")
    }

  } else {

    # Ensure all GRanges in GRangesList have correct widths
    assertthat::assert_that(
      all(sapply(tc_object,
                 function(gr) all(GenomicRanges::width(gr) == len_vec))),
      msg = paste("\n❌ All regions in each GRanges of tc_object (GRangesList) must have width", len_vec) # nolint: line_length_linter.
    )

    # Ensure the length of tc_object matches ctss_rse
    assertthat::assert_that(
      length(tc_object) == ncol(ctss_rse),
      msg = "\n❌ tc_object (GRangesList) must have the same length as ctss_rse" # nolint: line_length_linter.
    )

    # Ensure names match
    assertthat::assert_that(
      identical(names(tc_object), colnames(ctss_rse)),
      msg = "\n❌ tc_object (GRangesList) must have the same names as ctss_rse" # nolint: line_length_linter.
    )

    # Warn if any seqlevels differ
    ctss_sl <- GenomeInfoDb::seqlevels(ctss_rse)
    message("ctss_rse seqlevels: ", paste(ctss_sl, collapse = ", "))

    lapply(seq_along(tc_object), function(i) {
      tc_gr <- tc_object[[i]]
      tc_sl <- GenomeInfoDb::seqlevels(tc_gr)

      if (!setequal(tc_sl, ctss_sl)) {
        warning(sprintf("⚠️ seqlevels differ between sample '%s' and ctss_rse — pruning is NOT applied during validation.", # nolint: line_length_linter.
                        names(tc_object)[i]))
        message("  → tc_object seqlevels: ", paste(tc_sl, collapse = ", "))
      } else {
        message(sprintf("seqlevels match for sample '%s'",
                        names(tc_object)[i]))
      }
    })
  }
  message("Fnished validating tc_object.")
  return(TRUE)
}

#' Generate Sliding Windows for a Single Chromosome
#'
#' This function creates sliding windows over a single chromosome's
#' reduced tag cluster regions (as `GRanges`), typically generated from
#' extended CTSS clusters. Each window is centered at regular intervals
#' and then expanded on both ends by a specified distance.
#'
#' It is designed to be called internally for each chromosome by
#' `tc_sliding_window()` and should not require strand specificity.
#'
#' @param gr_per_chr A `GRanges` object containing ranges
#' for a single chromosome,
#' representing reduced and uniformly extended tag clusters (TCs).
#' @param sld_by Integer. Sliding window step size in base pairs (default: 20).
#' @param ext_dis Integer. Number of base pairs to expand the window
#' on both sides (default: 200).
#'
#' @return A `GRanges` object containing all expanded sliding windows
#' for the input chromosome.
#'
#' @examples
#' gr_chr <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = c(100, 200), end = c(150, 250))
#' )
#' windows <- tc_sliding_window_chr(gr_chr, sld_by = 20, ext_dis = 200)
#'
#' @import GenomicRanges
#' @import IRanges
#'
#' @keywords internal
tc_sliding_window_chr <- function(gr_per_chr,
                                  sld_by = 20,
                                  ext_dis = 200,
                                  log_file = log_file) {

  chr <- as.character(GenomicRanges::seqnames(gr_per_chr))[1]
  plc_log(sprintf("➡️ Starting sliding window for chromosome: %s",
                  chr),
          log_file, print_console = FALSE)

  # 1) Reduce overlaps within this chromosome
  collapsed_granges <- GenomicRanges::reduce(gr_per_chr)

  # 2) Generate sliding windows for all ranges in the collapsed GRanges
  sliding_granges_list <- lapply(seq_along(collapsed_granges),
                                 function(i) {
    start_pos <- GenomicRanges::start(collapsed_granges[i])
    end_pos <- GenomicRanges::end(collapsed_granges[i])

    # Adjust the sliding window start and end
    adjusted_start <- start_pos + ext_dis - sld_by
    adjusted_end <- end_pos - ext_dis + sld_by

    # Create a sequence of positions for sliding windows
    sliding_positions <- seq(adjusted_start, adjusted_end, by = sld_by)

    # Create a GRanges object for the sliding windows (width = 1)
    sliding_granges <- GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(collapsed_granges[i]),
      ranges = IRanges::IRanges(start = sliding_positions, width = 1),
      strand = GenomicRanges::strand(collapsed_granges[i])
    )

    # Store the original ranges (before expansion) as IRanges
    # in a "thick" metadata column
    sliding_granges$thick <- IRanges::IRanges(start = sliding_positions,
                                              width = 1)

    # Expand sliding_granges directly to save memory
    sliding_granges <- GenomicRanges::resize(sliding_granges,
                                             width = GenomicRanges::width(sliding_granges) + ext_dis, # nolint: line_length_linter.
                                             fix = "start")
    sliding_granges <- GenomicRanges::resize(sliding_granges,
                                             width = GenomicRanges::width(sliding_granges) + ext_dis, # nolint: line_length_linter.
                                             fix = "end")

    sliding_granges
  })

  # 3) Combine the list of GRanges objects into a single GRanges object
  sliding_granges <- base::do.call(c, sliding_granges_list)

  # 4) Remove out-of-bound sliding windows (start < 1)
  before <- length(sliding_granges)
  sliding_granges <- sliding_granges[start(sliding_granges) > 0]
  after <- length(sliding_granges)

  if (before > after) {
    plc_log(sprintf("⚠️ Removed %d out-of-bound windows (start < 1) in %s",
                    before - after, chr),
            log_file, print_console = FALSE)
  }

  plc_log(sprintf("✅ Finished sliding window for chromosome: %s",
                  chr),
          log_file, print_console = FALSE)

  rm(sliding_granges_list, collapsed_granges)
  gc()

  # Return the sliding GRanges object
  return(sliding_granges)
}

#' Perform Parallel Sliding Window Expansion on GRanges Tag Clusters
#'
#' @import GenomicRanges
#' @import IRanges
#' @importFrom future plan sequential multisession
#' @importFrom future.apply future_lapply
#' @importFrom assertthat assert_that
#' @importFrom parallel detectCores
#'
#' @export
tc_sliding_window <- function(granges_obj,
                              sld_by = 20,
                              ext_dis = 200,
                              log_file = log_file,
                              num_cores = NULL) {

  assertthat::assert_that(
    !is.null(granges_obj),
    length(granges_obj) > 0,
    msg = "\n❌ Error: Input granges_obj for tc_sliding_window() is NULL or empty. Cannot perform sliding window operation." # nolint: line_length_linter.
  )

  # 1) Ignore strand by setting all strands to "*"
  GenomicRanges::strand(granges_obj) <- "*"

  # 2) Split the GRanges object by chromosome (seqnames)
  gr_by_chr <- split(granges_obj, GenomicRanges::seqnames(granges_obj))
  num_jobs <- length(gr_by_chr)

  # 3) Determine the number of cores to use and global variable memory limit
  if (is.null(num_cores)) {
    num_cores <- max(1, min(25, parallel::detectCores() %/% 2))
  }
  num_workers <- min(num_cores, num_jobs)
  options(future.globals.maxSize = 4 * 1024^3)  # global variable memory 4GB limit # nolint: line_length_linter.
  future::plan(future::multisession, workers = num_workers)
  on.exit(future::plan(future::sequential))  # Reset future::plan() to default after execution # nolint: line_length_linter.

  plc_log(paste("Using", num_workers, "core(s) for parallel processing."),
          log_file, "INFO", print_console = FALSE)

  # 4) Run in parallel using future_lapply, passing required globals
  result_list <- future.apply::future_lapply(
    gr_by_chr,
    FUN = function(gr) {
      tc_sliding_window_chr(gr,
                            sld_by = sld_by,
                            ext_dis = ext_dis,
                            log_file = log_file)
    },
    future.seed = TRUE,
    future.globals = list(
      tc_sliding_window_chr = tc_sliding_window_chr,
      sld_by = sld_by,
      ext_dis = ext_dis,
      log_file = log_file
    )
  )

  # 5) Convert result into GRangesList and then unlist
  result_grl <- GenomicRanges::GRangesList(result_list)
  result_gr <- unlist(result_grl)

  # 6) Final sanity check
  assertthat::assert_that(
    !is.null(result_gr),
    length(result_gr) > 0, # nolint: line_length_linter.
    msg = "\n❌ Error: tc_sliding_window() returned NULL or an empty GRanges object." # nolint: line_length_linter.
  )

  return(result_gr)
}

#' Prepare Directory Structure for Profile Output
prep_profile_dir <- function(output_dir = ".",
                             profile_dir_name = "profile_output",
                             dir_name = c("metadata",
                                          "profiles",
                                          "profiles_subtnorm",
                                          "predictions")) {

  # Create the output dir and main dir
  new_path <- file.path(output_dir, "PRIMEloci_tmp", profile_dir_name)

  if (!file.exists(new_path)) {
    dir.create(new_path, recursive = TRUE, showWarnings = FALSE)

    # Create main directories and their subdirectories
    lapply(dir_name, function(main) {
      main_path <- file.path(new_path, main)
      dir.create(main_path)
    })

    message(paste0("New folder created:", new_path, "\n"))
  } else {
    message(paste0("Folder already exists:", new_path, "\n"))
  }
  return(new_path)
}

#' Convert SummarizedExperiment to GRanges with Assay Data
#'
#' This function converts a `SummarizedExperiment` object to a `GRanges` object
#' and adds a specified assay data column to the resulting `GRanges` object.
#'
#' @importFrom SummarizedExperiment rowRanges assay
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom magrittr %>%
cast_rse_to_granges <- function(rse,
                                assay = "counts",
                                coln_assay = 1,
                                colname = "score") {
  # Extract the row ranges from the SummarizedExperiment object
  gr <- SummarizedExperiment::rowRanges(rse) %>% GenomicRanges::GRanges() # nolint: pipe_operator_linter

  # Extract the assay data
  assay_data <- SummarizedExperiment::assay(rse, assay)

  # Assign the assay data to the specified column name in the GRanges object
  GenomicRanges::mcols(gr)[[colname]] <- assay_data[, coln_assay]

  return(gr)
}

#' Convert Strand Information to No Strand for GRanges Object
#'
#' @importFrom GenomicRanges strand
convert_strand_to_nostrand_gr <- function(region_gr) {
  GenomicRanges::strand(region_gr) <- "*"
  region_gr
}

#' Remove Metadata and Duplicate Genomic Ranges
#'
#' This function takes a `GRanges` object, removes all metadata,
#' and then eliminates duplicate genomic ranges based on sequence names,
#' start and end positions, and strand information.
#'
#' @importFrom GenomicRanges GRanges seqnames ranges strand start end duplicated
#' @importFrom IRanges IRanges
remove_metadata_and_duplicates <- function(gr) {
  # Remove metadata columns by creating a new GRanges object without metadata
  gr_no_metadata <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(gr),
    ranges = IRanges::IRanges(start = GenomicRanges::start(gr),
                              end = GenomicRanges::end(gr)),
    strand = GenomicRanges::strand(gr)
  )

  # Identify duplicated ranges based on seqnames, ranges, and strand
  duplicated_indices <- GenomicRanges::duplicated(gr_no_metadata)

  # Subset the GRanges object to keep only unique ranges
  unique_gr <- gr_no_metadata[!duplicated_indices]

  unique_gr
}

# Check presence and format of rownames in a profile matrix
check_valid_profile_rownames <- function(profile_matrix, chr_name, log_file) {
  rn <- rownames(profile_matrix)

  if (is.null(rn)) {
    msg <- sprintf("❌ No rownames found in heatmapData output for chromosome %s.", # nolint: line_length_linter.
                   chr_name)
    plc_log(msg, log_file, level = "❌ ERROR")
    stop(msg)
  }

  bad_rn <- rn[!grepl("^[^:]+:\\d+-\\d+;.$", rn)]

  if (length(bad_rn) > 0) {
    msg <- sprintf("❌ %d rownames have invalid format in heatmapData output for chromosome %s. Example: %s", # nolint: line_length_linter.
                   length(bad_rn), chr_name, bad_rn[1])
    plc_log(msg, log_file, level = "❌ ERROR")
    stop(msg)
  }
}

# Extract chr, start, end, and strand from row names like "chr:start-end;strand"
extract_rowname_components <- function(row_names) {
  str_parts <- strsplit(row_names, "[:]|-|;")

  # Convert to matrix (columns: chr, start, end, strand)
  components <- do.call(rbind, str_parts)

  if (ncol(components) != 4) {
    stop("❌ Invalid row name format. Expected 'chr:start-end;strand'")
  }

  components
}

# Create GRanges object from row names like "chr:start-end;strand"
create_granges_from_rownames <- function(row_names_str) {
  components <- extract_rowname_components(row_names_str)

  GenomicRanges::GRanges(
    seqnames = components[, 1],
    ranges = IRanges::IRanges(
      start = as.integer(components[, 2]),
      end   = as.integer(components[, 3])
    ),
    strand = components[, 4]
  )
}

# Convert strand suffix in row names to "*"
convert_rowname_to_nostrand <- function(strand_str) {
  gsub("[+-]$", "*", strand_str)
}

# Modify row names to no-strand if both strands match
modify_profile_rownames <- function(profiles, count_profiles) {
  plus_names <- convert_rowname_to_nostrand(rownames(count_profiles$`*`$`+`))
  minus_names <- convert_rowname_to_nostrand(rownames(count_profiles$`*`$`-`))

  if (identical(sort(plus_names), sort(minus_names))) {
    rownames(profiles) <- plus_names
  } else {
    warning("⚠️ Row names from + and - strands are not aligned after strand conversion.") # nolint: line_length_linter.
  }

  profiles
}

# Combine plus and minus strand profiles
combine_plus_minus_profiles <- function(count_profiles, len_vec) {
  plus <- count_profiles$`*`$`+`
  minus <- count_profiles$`*`$`-`

  combined <- cbind(plus, minus)

  colnames(combined) <- c(
    paste0("Plus_", seq_len(len_vec)),
    paste0("Minus_", seq_len(len_vec))
  )

  combined <- modify_profile_rownames(combined, count_profiles)

  combined
}

# Vectorized normalized strand subtraction on a matrix
strands_norm_subtraction <- function(mat, len_vec) {
  if (!inherits(mat, "matrix") && !inherits(mat, "Matrix")) {
    stop("❌ Input must be a dense matrix or sparse Matrix object.")
  }

  if (nrow(mat) == 0 || ncol(mat) < 2 * len_vec) {
    stop("❌ Matrix has insufficient dimensions for strand subtraction.")
  }

  # Extract plus and minus halves
  plus <- mat[, 1:len_vec, drop = FALSE]
  minus <- mat[, (len_vec + 1):(2 * len_vec), drop = FALSE]

  max_val <- sparseMatrixStats::rowMaxs(cbind(abs(plus), abs(minus)))
  norm_mat <- (plus - minus) / max_val

  norm_mat
}

# Apply normalized strand subtraction and assign position labels
strands_norm_subtraction_all <- function(windows,
                                         ext_dis,
                                         len_vec,
                                         output_sparse = TRUE) {
  if (is.null(dim(windows))) {
    stop("❌ `windows` must be a dense matrix or sparse matrix object.")
  }

  norm_mat <- strands_norm_subtraction(windows, len_vec)

  if (!output_sparse) {
    norm_mat <- as.matrix(norm_mat)
  }

  colnames(norm_mat) <- paste0("Pos", seq_len(ncol(norm_mat)) - ext_dis - 1)

  return(norm_mat)
}

# Save profile or metadata to file
# (csv, parquet, or npz with fallback and logging)
plc_save_to_file <- function(data,
                             suffix,
                             file_type,
                             chr_name,
                             file_path,
                             log_file) {
  full_file_path <- paste0(file_path, "_", chr_name, suffix)

  is_sparse <- inherits(data, "dgCMatrix") || inherits(data, "dgRMatrix")
  actual_file_type <- file_type

  if (file_type == "npz" && !is_sparse) {
    msg <- sprintf("⚠️ Requested npz output for %s, but data is not sparse. Saving as Parquet instead.", suffix) # nolint: line_length_linter.
    message(msg)
    plc_log(msg, log_file, level = "⚠️ WARN")
    actual_file_type <- "parquet"
  }

  options(scipen = 999)  # turn off scientific notation globally (during save)

  tryCatch({
    if (actual_file_type == "csv") {
      utils::write.csv(as.data.frame(data),
                       paste0(full_file_path, ".csv"),
                       row.names = FALSE)

    } else if (actual_file_type == "parquet") {
      arrow::write_parquet(as.data.frame(data),
                           paste0(full_file_path, ".parquet"))

    } else if (actual_file_type == "npz") {
      scipy <- reticulate::import("scipy.sparse", delay_load = TRUE)

      py_matrix <- if (inherits(data, "dgCMatrix")) {
        scipy$csc_matrix(reticulate::tuple(data@x, data@i, data@p),
                         shape = reticulate::tuple(data@Dim[1], data@Dim[2]))
      } else {
        scipy$csr_matrix(reticulate::tuple(data@x, data@j, data@p),
                         shape = reticulate::tuple(data@Dim[1], data@Dim[2]))
      }

      npz_path <- paste0(full_file_path, ".npz")
      scipy$save_npz(npz_path, py_matrix)

    } else {
      stop("❌ Unsupported file_type. Use 'csv', 'parquet', or 'npz'.")
    }

  }, error = function(e) {
    msg <- sprintf("❌ Failed to save %s for chromosome %s: %s",
                   suffix, chr_name, e$message)
    plc_log(msg, log_file, level = "❌ ERROR")
    stop(msg)
  })
}

#' Profile CTSS counts over strand-merged sliding windows
#' for a single chromosome.
#'
#' Internal function used by PRIMEloci to compute normalized
#' and raw signal profiles,
#' along with metadata including total signal (sum_count) per window.
#'
#' @param current_region_gr GRanges of sliding windows for a single chromosome.
#' @param filtered_ctss_gr GRanges of CTSS positions for the same chromosome.
#' @param chr_name Chromosome name.
#' @param file_path Output folder path.
#' @param output_file_prefix Output file prefix (excluding chr/suffix/ext).
#' @param ext_dis Extension distance used (to recover window length).
#' @param save_count_profiles Logical. Save raw count profiles as well?
#' @param file_type Output format: "csv", "parquet", or "npz".
#' @param log_file Log file path.
#'
#' @export
#'
#' @return A list with chromosome name and processing status.
plc_profile_chr <- function(current_region_gr,
                                  filtered_ctss_gr,
                                  chr_name,
                                  file_path,
                                  output_file_prefix,
                                  ext_dis,
                                  save_count_profiles = FALSE,
                                  file_type = "parquet",
                                  log_file) {

  plc_log(sprintf("➡️ Processing chromosome: %s", chr_name),
          log_file, print_console = FALSE)

  if (length(filtered_ctss_gr) == 0) {
    msg <- sprintf("⚠️ Skipping chromosome %s: No CTSS data found.",
                   chr_name)
    plc_log(msg, log_file)
    return(list(chr_name = chr_name, status = "Skipped: No CTSS"))
  }

  if (is.null(current_region_gr) || length(current_region_gr) == 0) {
    msg <- sprintf("⚠️ Skipping chromosome %s: No valid genomic regions found.",
                   chr_name)
    plc_log(msg, log_file)
    return(list(chr_name = chr_name, status = "Skipped: No regions"))
  }

  # Pre-processing
  current_region_gr <- convert_strand_to_nostrand_gr(current_region_gr)
  current_region_gr <- remove_metadata_and_duplicates(current_region_gr)

  # Generate sparse matrix profile
  count_profiles <- suppressMessages(
    heatmapData(current_region_gr,
                filtered_ctss_gr,
                sparse = TRUE)
  )

  check_valid_profile_rownames(count_profiles$`*`$`+`, chr_name, log_file)
  check_valid_profile_rownames(count_profiles$`*`$`-`, chr_name, log_file)

  len_vec <- ext_dis * 2 + 1

  # Combine plus/minus strand profiles
  combined_count_profiles <- combine_plus_minus_profiles(count_profiles,
                                                         len_vec)

  # Compute subtraction normalization
  output_sparse <- file_type == "npz"
  combined_subtnorm_profiles <- strands_norm_subtraction_all(combined_count_profiles, # nolint: line_length_linter.
                                                             ext_dis,
                                                             len_vec,
                                                             output_sparse = output_sparse) # nolint: line_length_linter.

  # check to confirm order of rownames
  if (!identical(rownames(combined_count_profiles),
                 rownames(combined_subtnorm_profiles))) {
    stop("❌ Rownames mismatch between count and subtnorm profiles.")
  }

  # Create metadata
  check_valid_profile_rownames(combined_count_profiles, chr_name, log_file)
  combined_count_metadata <- create_granges_from_rownames(rownames(combined_count_profiles)) # nolint: line_length_linter.
  if (length(combined_count_metadata) != nrow(combined_count_profiles)) {
    stop("❌ Length mismatch between matrix rows and generated GRanges.")
  }
  combined_count_metadata$sum_count <- Matrix::rowSums(combined_count_profiles)
  metadata_file_type <- if (file_type == "npz") "csv" else file_type

  # All saving

  # Save normalized profiles
  plc_save_to_file(combined_subtnorm_profiles,
                   "_profiles_subtnorm",
                   file_type,
                   chr_name,
                   file.path(file_path,
                             "profiles_subtnorm",
                             output_file_prefix),
                   log_file)
  rm(combined_subtnorm_profiles)

  if (save_count_profiles) {
    if (output_sparse) {
      plc_save_to_file(combined_count_profiles,
                       "_profiles",
                       file_type,
                       chr_name,
                       file.path(file_path,
                                 "profiles",
                                 output_file_prefix),
                       log_file)
    } else {
      plc_save_to_file(as.matrix(combined_count_profiles),
                       "_profiles",
                       file_type,
                       chr_name,
                       file.path(file_path,
                                 "profiles",
                                 output_file_prefix),
                       log_file)
    }
  }

  # Save metadata as parquet (if npz) or user-specified format
  plc_save_to_file(as.data.frame(combined_count_metadata),
                   "_metadata",
                   metadata_file_type,
                   chr_name,
                   file.path(file_path,
                             "metadata",
                             output_file_prefix),
                   log_file)

  # Cleanup
  rm(combined_count_metadata)
  rm(combined_count_profiles)
  gc()

  plc_log(sprintf("✅ Successfully processed chromosome: %s",
                  chr_name),
          log_file, print_console = FALSE)

  return(list(chr_name = chr_name, status = "Processed"))
}

#' Process profiles for each column in the CTSS dataset and save results
#'
#' This function processes the profiles for each column
#' in a RangedSummarizedExperiment object (`ctss_rse`),
#' computes the count profiles, and calls a chromosome-specific function
#' to handle individual chromosomes and save the results.
#'
#' @param ctss_rse A RangedSummarizedExperiment object containing CTSS counts.
#' @param regions_gr A GRanges or GRangesList object containing genomic regions.
#' @param output_dir The directory where output files will be saved.
#' @param profile_dir_name The directory name for the output.
#' @param ext_dis Integer value for extending the distance
#' used in profile computations.
#' @param python_path The path to the Python executable.
#' If NULL, the default Python environment will be used.
#' @param addtn_to_filename A string to add to the output filename.
#' @param save_count_profiles Logical flag indicating
#' whether count profiles should be saved.
#' @param file_type The output file format.
#' One of: "parquet" (default), "csv", or "npz".
#' @param num_cores The number of cores to use for parallel processing.
#' @param log_file_name The name of the log file
#' (default: "plc_log_profile").
#'
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom SummarizedExperiment colnames
#' @importFrom arrow write_parquet
#' @importFrom parallel detectCores
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom assertthat assert_that
#' @importFrom reticulate use_python
#' @export
plc_profile <- function(ctss_rse,
                              regions_gr,
                              output_dir,
                              profile_dir_name,
                              file_type = "npz",
                              python_path = NULL,
                              addtn_to_filename = "",
                              save_count_profiles = FALSE,
                              num_cores = NULL,
                              log_file = "./PRIMEloci_profile.log",
                              ext_dis) {

  # --- Assertions ---
  assertthat::assert_that(inherits(ctss_rse, "RangedSummarizedExperiment"),
                          msg = "❌ ctss_rse must be a RangedSummarizedExperiment object.") # nolint: line_length_linter.
  assertthat::assert_that(inherits(regions_gr, c("GRanges",
                                                 "GRangesList",
                                                 "CompressedGRangesList")),
                          msg = "❌ regions_gr must be a GRanges, GRangesList, or CompressedGRangesList object.") # nolint: line_length_linter.
  assertthat::assert_that(!is.null(regions_gr) && length(regions_gr) > 0,
                          msg = "❌ regions_gr must not be NULL or empty.")
  assertthat::assert_that(ext_dis %% 1 == 0,
                          msg = "❌ ext_dis must be an integer.")
  assertthat::assert_that(file_type %in% c("parquet", "csv", "npz"),
                          msg = "❌ file_type must be 'parquet', 'csv', or 'npz'.") # nolint: line_length_linter.

  # --- Setup ---
  file_path <- prep_profile_dir(output_dir = output_dir,
                                profile_dir_name = profile_dir_name)

  num_jobs <- length(GenomicRanges::seqnames(regions_gr))   
  if (is.null(num_cores)) {
    num_cores <- max(1, min(25, parallel::detectCores() %/% 2))
  }
  num_workers <- min(num_cores, num_jobs)
  options(future.globals.maxSize = 4 * 1024^3)  # global variable memory 4GB limit # nolint: line_length_linter.
  future::plan(future::multisession, workers = num_workers)
  on.exit(future::plan(future::sequential))  # Reset future::plan() to default after execution # nolint: line_length_linter.

  plc_log(paste("Using", num_workers, "core(s) for parallel processing."),
          log_file, "INFO", print_console = FALSE)

  # --- Sample Loop ---
  for (i in seq_along(SummarizedExperiment::colnames(ctss_rse))) {
    sample_name <- SummarizedExperiment::colnames(ctss_rse)[i]
    start_time <- Sys.time()

    # region_gr
    current_region_gr <- if (inherits(regions_gr,
                                      c("GRangesList",
                                        "CompressedGRangesList"))) {
      regions_gr[[i]]
    } else {
      regions_gr
    }

    # ctss_gr
    ctss_gr <- cast_rse_to_granges(ctss_rse, assay = "counts", coln_assay = i)
    regions_list <- split(current_region_gr,
                          GenomicRanges::seqnames(current_region_gr))

    # prepare output filename
    output_file_prefix <- paste0(sample_name, addtn_to_filename)

    results <- future.apply::future_lapply(names(regions_list),
                                           function(chr_name) {
      tryCatch({

        # 🔧 Set Python path for this child process
        reticulate::use_python(python_path, required = TRUE)

        # Filter ctss_gr for the current chromosome
        filtered_ctss_gr <- ctss_gr[as.character(GenomicRanges::seqnames(ctss_gr)) == chr_name] # nolint: line_length_linter.

        # Call the chromosome-specific function
        plc_profile_chr(
          regions_list[[chr_name]],
          filtered_ctss_gr,
          chr_name,
          file_path,
          output_file_prefix,
          ext_dis,
          save_count_profiles = save_count_profiles,
          file_type = file_type,
          log_file = log_file
        )
        list(chr_name = chr_name, status = "Processed")
      }, error = function(e) {
        plc_log(sprintf("Error processing chromosome %s for sample %s: %s",
                        chr_name, sample_name, e$message),
                log_file, level = "❌ ERROR")
        list(chr_name = chr_name, status = "Failed", error = e$message)
      })
    }, future.seed = TRUE)

    failed <- results[sapply(results, function(x) x$status == "Failed")]
    if (length(failed) > 0) {
      plc_log("The following chromosomes encountered errors:",
              log_file,
              level = "❌ ERROR")
      for (fail in failed) {
        plc_log(sprintf("Chromosome: %s | Error: %s",
                        fail$chr_name, fail$error),
                log_file,
                level = "❌ ERROR")
      }
    } else {
      plc_log("✅ All chromosomes processed successfully!", log_file)
    }

    plc_log(sprintf("✅ Finished processing sample: %s |  ✖️ Failed: %d/%d chromosomes", # nolint: line_length_linter.
                    sample_name, length(failed), length(regions_list)),
            log_file)

    runtime <- difftime(Sys.time(), start_time, units = "mins")
    plc_log(sprintf("⏱️ Time taken for sample %s: %.2f minutes",
                    sample_name, as.numeric(runtime)),
            log_file)

  }

  invisible(TRUE)

}

#' Load a BED file and validate its columns
#'
#' This function reads a BED file into a `data.table`
#' and checks that it contains the required columns:
#' 'chrom', 'chromStart', 'chromEnd', 'strand', and 'score'.
#'
#' @param input_bed Character. The path to the input BED file.
#'
#' @return A `data.table` containing the BED file data.
#'
#' @importFrom data.table fread
#' @importFrom assertthat assert_that
#' @export
load_bed_file <- function(input_bed) {
  bed_file <- read.table(input_bed,
                         header = TRUE,
                         sep = "\t",
                         stringsAsFactors = FALSE)
  required_cols <- c("chrom", "chromStart", "chromEnd", "strand", "score")
  assertthat::assert_that(all(required_cols %in% colnames(bed_file)),
                          msg = "The BED file must contain 'chrom', 'chromStart', 'chromEnd', 'strand', and 'score' columns.") # nolint: line_length_linter.
  bed_file
}

#' Create a GRanges object from a BED data.table
#'
#' This function converts a `data.table`
#' containing BED file data into a `GRanges` object.
#' The required columns are extracted and
#' used to define the `GRanges` object, and the remaining
#' columns are added as metadata.
#'
#' @param bed_file A `data.table` containing the BED file data.
#'
#' @return A `GRanges` object.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
create_granges_from_bed <- function(bed_file) {
  gr <- GenomicRanges::GRanges(seqnames = bed_file$chrom,
                               ranges = IRanges::IRanges(start = bed_file$chromStart + 1, # nolint: line_length_linter.
                                                         end = bed_file$chromEnd), # nolint: line_length_linter.
                               strand = bed_file$strand)
  S4Vectors::mcols(gr) <- bed_file[, !(names(bed_file) %in% c("chrom",
                                                              "chromStart",
                                                              "chromEnd",
                                                              "strand"))]
  gr
}

#' Selectively merge overlapping cores based on score difference.
#'
#' This function takes a GRanges object of genomic cores, sorts them
#' by descending score, and merges overlapping cores if their score
#' difference is within a specified threshold.
#'
#' @param core_gr A GRanges object containing genomic core ranges.
#' @param score_diff A numeric value specifying the maximum allowed
#' score difference for merging overlapping cores.
#' @return A GRanges object with merged core regions and added
#' metadata, including the thick position and maximum score.
#' @importFrom GenomicRanges GRanges reduce findOverlaps
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom S4Vectors mcols subjectHits
selective_merge_cores <- function(core_gr, score_diff) {
  # Sort cores by descending score
  core_gr <- core_gr[order(-core_gr$score)]

  # Initialize final GRanges and metadata containers
  merged_cores <- GenomicRanges::GRanges()
  thick_vals <- IRanges::IRanges()
  max_scores <- numeric(0)

  while (length(core_gr) > 0) {
    # Take the highest-ranked core
    x <- core_gr[1]
    score_x <- x$score
    thick_x <- x$thick

    ## 1. Identify overlapping cores to the top core
    overlaps <- GenomicRanges::findOverlaps(x, core_gr)
    overlap_set <- core_gr[S4Vectors::subjectHits(overlaps)]

    ## 2. Deviate at most "score_diff aka d" in score from the top core
    merge_candidates <- overlap_set[overlap_set$score >= score_x - score_diff]

    if (length(merge_candidates) > 1) {
      # 3. Merge the merge_candidates
      merged_region <- GenomicRanges::reduce(merge_candidates)
      thick_vals <- c(thick_vals, thick_x)
      max_scores <- c(max_scores, score_x)
      merged_cores <- c(merged_cores, merged_region)

      # 4. Remove from the set of cores,
      # those that overlap with the merged cores
      core_gr <- IRanges::subsetByOverlaps(core_gr, merged_cores, invert = TRUE)

    } else {
      # 5. If there is no merge, keep the top core
      # and remove overlaps from the set
      thick_vals <- c(thick_vals, thick_x)
      max_scores <- c(max_scores, score_x)
      merged_cores <- c(merged_cores, x)

      core_gr <- core_gr[-S4Vectors::subjectHits(overlaps)]
    }


  }

  # Attach metadata
  mcols(merged_cores)$thick <- thick_vals
  mcols(merged_cores)$max_score <- max_scores

  merged_cores
}

#' Extract sample label from input_basename (internal)
#'
#' Returns the substring between 'pred_all_' and '_combined' if both exist;
#' otherwise returns 'PRIMEloci'.
extract_sample_label <- function(input_basename) {
  pattern <- "pred_all_(.*?)_combined"
  matches <- stringr::str_match(input_basename, pattern)

  if (!is.na(matches[2])) {
    paste0("PRIMEloci_", matches[2])
  } else {
    "PRIMEloci"
  }
}

#' Write a GRanges object to a BED file for coreovlwith-d.
#'
#' This function converts a GRanges object to a data frame
#' and writes it to a BED file, including metadata such as
#' thick position and maximum score.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @importFrom data.table fwrite
write_granges_to_bed_coreovlwithd <- function(gr,
                                              output_dir,
                                              input_basename,
                                              output_basename) {
  output_bed <- file.path(output_dir, paste0(output_basename, ".bed"))

  bed_df <- as.data.frame(gr)

  if (!"thick" %in% colnames(mcols(gr)) || !"max_score" %in% colnames(mcols(gr))) { # nolint: line_length_linter.
    stop("'thick' or 'max_score' column not found in the metadata.")
  }

  bed_df$thick <- as.numeric(start(gr$thick))
  bed_df$start <- bed_df$start - 1  # Convert start to 0-based for BED format
  bed_df <- within(bed_df, {
    thickStart <- thick - 1
    thickEnd <- thick
    name <- paste0(input_basename, "-", seq_len(nrow(bed_df)))
  })

  # Reorder and rename columns
  bed_df <- bed_df[, c("seqnames", "start", "end",
                       "name", "max_score", "strand",
                       "thickStart", "thickEnd")]
  data.table::setnames(bed_df,
                       c("seqnames", "start", "end", "strand"),
                       c("chrom", "chromStart", "chromEnd", "strand"))
  # Write to BED file
  data.table::fwrite(bed_df,
                     file = output_bed,
                     sep = "\t",
                     quote = FALSE,
                     col.names = FALSE)
  cat("Reduced GRanges object saved to", output_bed, "\n")
}

#' Collapse core regions from a BED file with score-based filtering
#'
#' This function reads a BED file, filters regions by score,
#' resizes them to a fixed core width, and merges overlapping cores selectively
#' based on score differences. It supports chromosome-wise parallel processing
#' and optional output to a BED file with logging throughout.
#'
#' @param bed_file Path to the input BED file.
#' @param score_threshold Numeric threshold for filtering regions
#' by score from 0-1. Default is 0.75.
#' @param score_diff Maximum allowed score difference
#' between overlapping regions for merging. Default is 0.1.
#' @param core_width Width to which each region should be resized (centered).
#' Default is 151.
#' @param return_gr Logical. If TRUE,
#' returns the final collapsed GRanges object. Default is TRUE.
#' @param output_dir Optional directory path to write BED output.
#' If NULL or FALSE, no output is written.
#' @param num_cores Number of CPU cores to use.
#' If NULL, will use half of available cores (up to 25). Default is NULL.
#' @param log_file Path to the log file. Default is "coreovl_with_d.log".
#'
#' @return If `return_gr = TRUE`, returns a `GRanges` object
#' containing the collapsed regions. Otherwise, returns `NULL`.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Loads and filters BED entries by score.
#'   \item Resizes regions to a fixed core width centered on each region.
#'   \item Merges overlapping cores selectively per chromosome
#'         using a parallel backend.
#'   \item Logs progress and optionally writes output to BED format.
#' }
#'
#' If no regions pass the score filter, the function exits early with a warning.
#'
#' @import GenomicRanges
#' @import S4Vectors
#' @import IRanges
#' @import assertthat
#' @importFrom parallel detectCores
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom stringr str_replace_all
#' @importFrom tools file_path_sans_ext
#' @importFrom magrittr %>%
coreovl_with_d <- function(bed_file,
                           score_threshold = 0.75,
                           score_diff = 0.1,
                           core_width = 151,
                           return_gr = TRUE,
                           output_dir = NULL,
                           num_cores = NULL,
                           log_file = "coreovl_with_d.log") {

  # Initialize logging
  plc_log("Starting coreovl_with_d()", log_file)

  assert_that(file.exists(bed_file),
              msg = "❌ Input BED file not found.")

  assert_that(is.number(score_threshold),
              score_threshold >= 0,
              score_threshold < 1,
              msg = "❌ score_threshold must be a number between 0 and 1.")

  assert_that(is.number(score_diff),
              score_diff >= 0,
              score_diff < score_threshold,
              msg = "❌ score_diff must be non-negative and less than score_threshold.") # nolint: line_length_linter.

  assert_that(is.count(core_width),
              msg = "❌ core_width must be a positive integer.")

  if (!is.null(output_dir)) {
    assert_that(is.string(output_dir),
                msg = "❌ output_dir must be a character string.")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      message(sprintf("✅ Created output directory: %s", output_dir))
    }
  }

  if (!is.null(num_cores)) {
    assert_that(is.count(num_cores),
                msg = "❌ num_cores must be a positive integer.")
  }

  # Load and prepare data
  bed <- load_bed_file(bed_file)
  gr <- create_granges_from_bed(bed)

  # Filter by score threshold
  filtered_gr <- gr[gr$score >= score_threshold]

  if (length(filtered_gr) == 0) {
    plc_log("No entries passed the score threshold. Exiting.",
            log_file, level = "⚠️ WARN")
    return(NULL)
  }

  chr_list <- unique(as.character(GenomicRanges::seqnames(filtered_gr)))
  num_jobs <- length(chr_list)
  error_messages <- list()  # store errors per chromosome

  # Setup parallelism

  if (is.null(num_cores)) {
    num_cores <- max(1, min(25, parallel::detectCores() %/% 2))
  }
  num_workers <- min(num_cores, num_jobs)
  options(future.globals.maxSize = 4 * 1024^3)  # global variable memory 4GB limit # nolint: line_length_linter.
  future::plan(future::multisession, workers = num_workers)
  on.exit(future::plan(future::sequential))  # Reset future::plan() to default after execution # nolint: line_length_linter.

  plc_log(paste("Using", num_workers, "core(s) for parallel processing."),
          log_file, "INFO", print_console = FALSE)

  # Process each chromosome in parallel

  collapsed_gr_list <- future.apply::future_lapply(chr_list, function(chr) {
    tryCatch({

      chr_vec <- as.character(GenomicRanges::seqnames(filtered_gr))
      chr_gr <- filtered_gr[chr_vec == chr]

      core_gr <- GenomicRanges::resize(chr_gr,
                                       width = core_width,
                                       fix = "center")
      core_gr$thick <- GenomicRanges::start(core_gr) + floor(core_width / 2)
      result <- selective_merge_cores(core_gr, score_diff)

      plc_log(sprintf("✅ Finished processing chromosome %s", chr),
              log_file, print_console = FALSE)

      result
    }, error = function(e) {
      error_messages[[chr]] <<- conditionMessage(e)
      plc_log(
        sprintf("❌ Error processing chromosome %s: %s",
                chr, conditionMessage(e)),
        log_file, level = "❌ ERROR", print_console = FALSE
      )
      NULL
    })
  })

  if (length(error_messages) > 0) {
    failed_chr <- names(error_messages)
    plc_log(
      paste("⚠️ Skipped", length(failed_chr),
            "chromosomes due to errors:",
            paste(failed_chr, collapse = ", ")),
      log_file,
      level = "⚠️ WARN"
    )
  }

  # Check again — do we have any GRanges?
  if (length(collapsed_gr_list) == 0) {
    plc_log("⚠️ No valid GRanges to collapse — all chromosomes failed or were empty.", # nolint: line_length_linter.
            log_file, level = "⚠️ WARN", print_console = TRUE)
    return(NULL)
  }

  # Collapse and validate type
  collapsed_gr <- do.call(c, collapsed_gr_list)

  if (!inherits(collapsed_gr, "GRanges")) {
    plc_log("❌ Combined result is not a valid GRanges object. Skipping sort.",
            log_file, level = "❌ ERROR", print_console = TRUE)
    return(NULL)
  }

  mcols(collapsed_gr) <- mcols(collapsed_gr)[, c("thick", "max_score"), drop = FALSE] # nolint: line_length_linter.

  collapsed_gr <- GenomeInfoDb::sortSeqlevels(collapsed_gr)
  collapsed_gr <- GenomicRanges::sort(collapsed_gr)

  # Output writing
  plc_log("Processing the output...", log_file)
  if (!is.null(output_dir)) {
    input_basename <- tools::file_path_sans_ext(basename(bed_file)) %>%
      stringr::str_replace_all("[^[:alnum:]]", "_")
    sample_label <- extract_sample_label(input_basename)

    # Append threshold and d to filename
    output_basename <- sprintf("%s_thresh%s_d%s",
                               input_basename,
                               format(score_threshold,
                                      digits = 2,
                                      scientific = FALSE),
                               format(score_diff,
                                      digits = 2,
                                      scientific = FALSE))

    plc_log(sprintf("Writing BED output: %s.bed", output_basename), log_file)

    write_granges_to_bed_coreovlwithd(
      collapsed_gr,
      output_dir,
      sample_label,
      output_basename
    )

    plc_log("✅ coreovl_with_d() finished successfully.", log_file)
  } else {
    plc_log("⚠️ skipping file writing at the end of coreovl_with_d step.",
            log_file, level = "⚠️ WARN", print_console = TRUE)
  }

  if (return_gr) return(collapsed_gr)

}

#' Find .bed files matching a partial name in a directory and log if none found
#'
#' @param dir Path to the directory to search.
#' @param partial_name Partial file name to match (not case-sensitive).
#' @param log_file Path to the log file (optional).
#' If provided, logs will be written.
#'
#' @return A character vector of matching .bed file paths (may be empty).
#' @export
find_bed_files_by_partial_name <- function(dir,
                                           partial_name,
                                           log_file = NULL) {
  pattern <- paste0("(?i)", partial_name, ".*\\.bed$")
  files <- list.files(
    path = dir,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(files) == 0) {
    msg <- sprintf("⚠️ No .bed files found in '%s' matching '%s'",
                   dir, partial_name)
    if (!is.null(log_file)) {
      plc_log(msg, log_file, level = "⚠️ WARN")
    } else {
      warning(msg)
    }
  }

  return(files)
}