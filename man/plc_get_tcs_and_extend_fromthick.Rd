% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PRIMEloci_utility.R
\name{plc_get_tcs_and_extend_fromthick}
\alias{plc_get_tcs_and_extend_fromthick}
\title{Get tag clusters and extend from thick positions}
\usage{
plc_get_tcs_and_extend_fromthick(ctss_rse, ext_dis = 200)
}
\arguments{
\item{ctss_rse}{A RangedSummarizedExperiment object containing CAGE data.}

\item{ext_dis}{An integer specifying the distance
to extend around thick positions (default is 200).}
}
\value{
A GenomicRanges::GRangesList object where each element corresponds to
        a tag cluster extended by \code{ext_dis} around the thick positions.
}
\description{
This function calculates tag clusters from
a RangedSummarizedExperiment object \code{ctss_rse}
and extends them by a specified distance \code{ext_dis}
around the thick positions.
}
\details{
The function iterates over columns of \code{ctss_rse},
calculates pooled counts,
subsets clusters based on a score threshold (> 0),
and clusters them unidirectionally.
It then extends each cluster by \code{ext_dis} base pairs
around the thick positions
and stores the results in a GRangesList object.
}
\examples{
# Example usage with a RangedSummarizedExperiment object
# ctss_rse <- ...  # Load or create your RangedSummarizedExperiment object
# ext_dis <- 200   # Define your extension distance
# result <- get_tagclusters_and_extend_fromthick(ctss_rse, ext_dis)

}
