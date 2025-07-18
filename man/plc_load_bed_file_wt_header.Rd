% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PRIMEloci_utility.R
\name{plc_load_bed_file_wt_header}
\alias{plc_load_bed_file_wt_header}
\title{Load a BED file and validate required columns}
\usage{
plc_load_bed_file_wt_header(input_bed, sum_count = FALSE)
}
\arguments{
\item{input_bed}{Character. Path to the input BED file.}

\item{sum_count}{Logical. Whether to require the 'sum_count' column
(default: FALSE).}
}
\value{
A `data.table` containing the BED file data with validated columns.
}
\description{
This function reads a BED file into a `data.table` and verifies that
it includes the necessary columns for downstream processing.
By default, it checks for the presence of 'chrom', 'chromStart',
'chromEnd', 'strand', and 'score'. If `sum_count = TRUE`, it also
requires the 'sum_count' column.
}
