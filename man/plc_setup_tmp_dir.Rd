% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PRIMEloci_utility.R
\name{plc_setup_tmp_dir}
\alias{plc_setup_tmp_dir}
\title{Set up temporary directory under output directory}
\usage{
plc_setup_tmp_dir(output_dir)
}
\arguments{
\item{output_dir}{A character string specifying the path
to the output directory.}
}
\value{
A character string with the full path to the temporary directory.
}
\description{
Creates a `PRIMEloci_tmp` folder under the given output directory
if it does not already exist.
}
