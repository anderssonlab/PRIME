#' Log messages to a file
#'
#' This function writes log messages to a specified log file,
#' appending a timestamp to each entry.
#'
#' @param message A character string containing the log message to write.
#' @param log_file A character string specifying
#' the file path where the log message should be written.
#'
#' @importFrom utils write.table
#'
profile_log_message <- function(message, log_file) {
  write.table(paste(Sys.time(), message),
              file = log_file,
              append = TRUE,
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
}



#' Process profiles for a single chromosome and save results
#'
#' This function processes the profiles for a specific chromosome
#' (GRanges object),
#' computes the count profiles, combines the results, and saves them in the
#' specified format (CSV or Parquet). The filename is provided
#' by the main function,
#' and this function appends the chromosome name to the filename.
#'
#' @param current_region_gr GRanges object for the current chromosome.
#' @param chr_name The name of the chromosome being processed.
#' @param filtered_ctss_gr A GRanges object containing
#' CTSS data filtered for the current chromosome.
#' @param file_path The base file path where the output will be saved.
#' @param base_file_name The base filename for the output files.
#' @param ext_dis Integer value for extending the distance
#' used in profile computations.
#' @param save_count_profiles Logical flag indicating
#' whether count profiles should be saved.
#' @param file_type The format in which to save the files ("csv" or "parquet").
#' @param log_file Path to the log file
#' where processing messages will be written.
#'
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom data.table fwrite
#' @importFrom arrow write_parquet
#' @importFrom PRIME heatmapData
#' @importFrom utils write.table
#'
PRIMEloci_profile_chr_2 <- function(current_region_gr,
                                    filtered_ctss_gr,
                                    chr_name,
                                    file_path,
                                    base_file_name,
                                    ext_dis,
                                    save_count_profiles,
                                    file_type,
                                    log_file) {
  profile_log_message(sprintf("🟢 Processing chromosome: %s",
                              chr_name),
                      log_file)

  # Ensure CTSS data is present
  if (length(filtered_ctss_gr) == 0) {
    profile_log_message(sprintf("⚠️ Skipping chromosome %s: No CTSS data available.", # nolint: line_length_linter.
                                chr_name),
                        log_file)
    return(list(chr_name = chr_name, status = "Skipped: No CTSS data"))
  }

  # Ensure genomic regions are valid before processing
  if (is.null(current_region_gr) || length(current_region_gr) == 0) {
    profile_log_message(sprintf("⚠️ Skipping chromosome %s: No valid genomic regions found.", # nolint: line_length_linter.
                                chr_name),
                        log_file)
    return(list(chr_name = chr_name, status = "Skipped: No valid regions"))
  }

  # Pre-processing: Convert strand and remove duplicates
  current_region_gr <- convert_strand_to_nostrand_gr(current_region_gr)
  current_region_gr <- remove_metadata_and_duplicates(current_region_gr)

  # Compute count profiles
  count_profiles <- heatmapData(current_region_gr, filtered_ctss_gr))
  rm(current_region_gr, filtered_ctss_gr)

  len_vec <- ext_dis * 2 + 1

  # Combine plus/minus strand profiles
  combined_count_profiles <- combine_plus_minus_profiles(count_profiles,
                                                         len_vec)
  rm(count_profiles)

  # Compute subtraction normalization
  combined_subtnorm_profiles <- strands_norm_subtraction_all(combined_count_profiles, # nolint: line_length_linter.
                                                             ext_dis,
                                                             len_vec)

  # Create metadata
  combined_count_metadata <- create_granges_from_rownames(rownames(combined_count_profiles)) # nolint: line_length_linter.
  sum_count <- data.frame(rowSums(combined_count_profiles))
  colnames(sum_count) <- "sum_count"
  combined_count_metadata$sum_count <- sum_count

  # Save profiles
  combined_subtnorm_profiles$rownames <- rownames(combined_subtnorm_profiles)
  save_to_file(combined_subtnorm_profiles, "_profiles_subtnorm",
               file_type, chr_name,
               file.path(file_path, "profiles_subtnorm", base_file_name))
  rm(combined_subtnorm_profiles)

  # Save count profiles if requested
  if (save_count_profiles) {
    combined_count_profiles$rownames <- rownames(combined_count_profiles)
    save_to_file(combined_count_profiles, "_profiles",
                 file_type, chr_name,
                 file.path(file_path, "profiles", base_file_name))
  }

  # Save metadata per chromosome
  combined_count_metadata$rownames <- rownames(combined_count_metadata)
  save_to_file(combined_count_metadata, "_metadata",
               file_type, chr_name,
               file.path(file_path, "metadata", base_file_name))
  rm(combined_count_metadata)
  rm(combined_count_profiles)

  # Garbage collection to free memory
  gc()

  profile_log_message(sprintf("✅ Successfully processed chromosome: %s",
                              chr_name),
                      log_file)
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
#' @param output_dir_name The directory name for the output.
#' @param ext_dis Integer value for extending the distance
#' used in profile computations.
#' @param addtn_to_filename A string to add to the output filename.
#' @param save_count_profiles Logical flag indicating
#' whether count profiles should be saved.
#' @param file_type The format in which to save the files ("csv" or "parquet").
#' @param num_cores The number of cores to use for parallel processing.
#' If not provided, it defaults to using all available cores minus one.
#'
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom SummarizedExperiment colnames
#' @importFrom arrow write_parquet
#' @importFrom parallel detectCores
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom utils write.table
#' @export
PRIMEloci_profile_2 <- function(ctss_rse,
                                regions_gr,
                                output_dir,
                                output_dir_name,
                                ext_dis,
                                addtn_to_filename = "",
                                save_count_profiles = FALSE,
                                file_type = "parquet",
                                num_cores = NULL) {

  # Prepare the output directory
  prep_profile_dir(output_dir = output_dir, output_dir_name = output_dir_name)

  # Detect available cores, leaving one free
  if (is.null(num_cores)) {
    num_cores <- max(1, parallel::detectCores() - 1)
  }

  # Set up parallel execution explicitly using future
  future::plan(future::multisession, workers = num_cores)

  # Create log directory and log file
  log_dir <- file.path(output_dir, "logs")
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  log_file <- file.path(log_dir, "PRIMEloci_log.txt")

  # Initialize log file
  write(sprintf("🔹 PRIMEloci Processing Started at %s\n",
                Sys.time()),
        file = log_file)

  print("PRIMEloci profile ...")

  # Process each column in ctss_rse
  for (i in seq_along(SummarizedExperiment::colnames(ctss_rse))) {
    sample_name <- SummarizedExperiment::colnames(ctss_rse)[i]
    profile_log_message(sprintf("🟢 Processing sample: %s",
                                sample_name),
                        log_file)
    print(sprintf("Processing sample: %s", sample_name))

    # Start time logging
    start_time <- Sys.time()
    profile_log_message(sprintf("⏳ Start time: %s", start_time), log_file)

    # Handle GRangesList and GRanges objects
    if (inherits(regions_gr, "GRangesList")) {
      current_region_gr <- regions_gr[[i]]
    } else if (inherits(regions_gr, "GRanges")) {
      current_region_gr <- regions_gr
    } else {
      stop("❌ Error: regions_gr is neither GRanges nor GRangesList.")
    }

    # Extract the GRanges for the current column of ctss_rse
    ctss_gr <- cast_rse_to_granges(ctss_rse, assay = "counts", coln_assay = i)

    # Split regions by chromosome
    regions_list <- split(current_region_gr,
                          GenomicRanges::seqnames(current_region_gr))

    # Generate the base file path
    file_path <- file.path(output_dir, output_dir_name)
    base_file_path <- paste0(sample_name, addtn_to_filename)

    # Parallel execution for each chromosome
    results <- future.apply::future_lapply(names(regions_list),
                                           function(chr_name) {
      tryCatch({
        # Filter ctss_gr for the current chromosome
        filtered_ctss_gr <- ctss_gr[GenomicRanges::seqnames(ctss_gr) == chr_name] # nolint: line_length_linter.

        # Call the chromosome-specific function
        PRIMEloci_profile_chr_2(
          current_region_gr = regions_list[[chr_name]],
          chr_name = chr_name,
          filtered_ctss_gr = filtered_ctss_gr,
          file_path = file_path,
          base_file_name = base_file_path,
          ext_dis = ext_dis,
          save_count_profiles = save_count_profiles,
          file_type = file_type,
          log_file = log_file
        )

        # If successful, return a success message
        return(list(chr_name = chr_name, status = "Processed"))

      }, error = function(e) {
        profile_log_message(sprintf("❌ Error processing chromosome %s for sample %s: %s", # nolint: line_length_linter.
                                    chr_name,
                                    sample_name,
                                    e$message),
                            log_file)
        return(list(chr_name = chr_name, status = "Failed", error = e$message))
      })
    }, future.seed = TRUE)

  # Summarize results after parallel execution
  failed_chromosomes <- results[sapply(results,
                                       function(x) x$status == "Failed")]

  # Log failures if any chromosomes failed
  if (length(failed_chromosomes) > 0) {
    profile_log_message("❌ The following chromosomes encountered errors:",
                        log_file)
    for (fail in failed_chromosomes) {
      profile_log_message(sprintf("Chromosome: %s | Error: %s",
                                  fail$chr_name,
                                  fail$error),
                          log_file)
    }
  } else {
    profile_log_message("✅ All chromosomes processed successfully!", log_file)
  }
  }
}

