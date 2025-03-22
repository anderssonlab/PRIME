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
                    print_console = TRUE) {
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
                level = "❌ERROR",
                print_console = print_console)
        NULL  # return NULL on error
      }
    ),
    warning = function(w) {
      plc_log(conditionMessage(w),
              log_file,
              level = "⚠️WARN",
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
      message("✅ seqlevels match between tc_object and ctss_rse")
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
        message(sprintf("✅ seqlevels match for sample '%s'",
                        names(tc_object)[i]))
      }
    })
  }
  message("✅ TC object validation passed successfully.")
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

    return(sliding_granges)
  })

  # 3) Combine the list of GRanges objects into a single GRanges object
  sliding_granges <- base::do.call(c, sliding_granges_list)

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

  # 3) Determine the number of cores to use and global variable memory limit
  if (is.null(num_cores)) {
    num_cores <- min(25, parallel::detectCores() %/% 2)
  }
  options(future.globals.maxSize = 4 * 1024^3)  # global variable memory 4GB limit # nolint: line_length_linter.
  future::plan(future::multisession, workers = num_cores)
  on.exit(future::plan(future::sequential))  # Reset future::plan() to default after execution # nolint: line_length_linter.

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
