% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PRIMEloci_utility.R
\name{plc_profile_chr}
\alias{plc_profile_chr}
\title{Profile CTSS counts over strand-merged sliding windows
for a single chromosome.}
\usage{
plc_profile_chr(
  current_region_gr,
  filtered_ctss_gr,
  chr_name,
  file_path,
  output_file_prefix,
  ext_dis,
  save_count_profiles = FALSE,
  file_type = "parquet"
)
}
\arguments{
\item{current_region_gr}{GRanges of sliding windows for a single chromosome.}

\item{filtered_ctss_gr}{GRanges of CTSS positions for the same chromosome.}

\item{chr_name}{Chromosome name.}

\item{file_path}{Output folder path.}

\item{output_file_prefix}{Output file prefix (excluding chr/suffix/ext).}

\item{ext_dis}{Extension distance used (to recover window length).}

\item{save_count_profiles}{Logical. Save raw count profiles as well?}

\item{file_type}{Output format: "csv", "parquet", or "npz".}
}
\value{
A list with chromosome name and processing status.
}
\description{
Internal function used by PRIMEloci to compute normalized
and raw signal profiles,
along with metadata including total signal (sum_count) per window.
}
