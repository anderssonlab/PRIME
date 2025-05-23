% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/noise.R
\name{estimateDivergentNoise}
\alias{estimateDivergentNoise}
\title{Estimate divergent loci noise in a given CTSS dataset by sampling and
quantifying regions of a fixed size across the genome.}
\usage{
estimateDivergentNoise(
  object,
  mask,
  mappable_minus,
  mappable_plus,
  map_frac = 0.5,
  win_size = 200,
  num_win = 1e+06,
  inputAssay = "counts",
  quantiles = c(0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999)
)
}
\arguments{
\item{object}{A \code{RangedSummarizedExperiment} object.}

\item{mask}{A \code{GRanges} object with the masked regions.}

\item{mappable_minus}{A \code{GRanges} object with the mappable regions for 
the minus strand.}

\item{mappable_plus}{A \code{GRanges} object with the mappable regions for 
the plus strand.}

\item{map_frac}{The fraction of bases that need to be mappable.}

\item{win_size}{The size of the windows to use.}

\item{num_win}{The number of windows to sample.}

\item{inputAssay}{The assay to use.}

\item{quantiles}{The quantiles to calculate.}
}
\value{
A matrix with the quantiles for each sample.
}
\description{
Estimate divergent loci noise in a given CTSS dataset by sampling and
quantifying regions of a fixed size across the genome.
}
