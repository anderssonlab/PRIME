#' Calculate the gene and CTSS complexity of a CTSS dataset.
#'
#' This function calculates the CTSS and gene complexity of a given CTSS
#' dataset.
#'
#' @param object A \code{RangedSummarizedExperiment} resulting from
#' CAGEfightR::quantifyCTSSs().
#' @param txModels A \code{TxDb} object containing the gene models.
#' @param step The step size for the subsampling.
#' @param CTSSunexpressed The count threshold for considering a CTSS to be
#' unexpressed.
#' @param geneunexpressed The count threshold for considering a gene to be
#' unexpressed.
#' @param minCTSSsupport The minimum number of samples a CTSS must be
#' expressed in to be considered.
#'
#' @return A data frame containing the number of CTSSs (column "numberCTSSs")
#' and genes (column "numberGenes") for each sample and target (column
#' "totalTags").
#'
#' @export
#'
#' @import CAGEfightR
#' @importFrom assertthat assert_that is.flag
#' @importFrom methods is
#'
#'
calcComplexity <- function(object, txModels, step = 1e6, CTSSunexpressed = 1,
                           geneunexpressed = 9, minCTSSsupport = 2) {
  assertthat::assert_that(
    methods::is(object, "RangedSummarizedExperiment"),
    methods::is(txModels, "TxDb"),
    is.numeric(step), length(step) == 1, step > 0,
    is.numeric(CTSSunexpressed), length(CTSSunexpressed) == 1, CTSSunexpressed >= 0,
    is.numeric(geneunexpressed), length(geneunexpressed) == 1, geneunexpressed >= 0,
    is.numeric(minCTSSsupport), length(minCTSSsupport) == 1, minCTSSsupport >= 1
  )

  object <- suppressWarnings(calcTotalTags(object,
    inputAssay = "counts"
  ))
  object <- suppressMessages(suppressWarnings(
    assignGeneID(object,
      geneModels = txModels,
      outputColumn = "geneID"
    )
  ))

  targets <- seq(step, max(object$totalTags), by = step)

  message("Subsampling counts")
  res <- c(
    list(data.frame(
      target = 0, sample = colnames(object), totalTags = 0,
      numberCTSSs = 0, numberGenes = 0
    )),
    lapply(targets, function(t) {
      message(t)

      x <- subsampleTarget(object, "counts", t)
      x <- suppressWarnings(calcTotalTags(x, inputAssay = "counts"))
      if (minCTSSsupport > 1) {
        x <- suppressMessages(subsetBySupport(x,
          unexpressed = CTSSunexpressed,
          minSamples = minCTSSsupport
        ))
      }

      x <- suppressWarnings(calcNumberCTSSs(x,
        inputAssay = "counts",
        unexpressed = CTSSunexpressed
      ))
      x <- suppressWarnings(calcNumberGenes(x,
        txModels,
        inputAssay = "counts",
        unexpressed = geneunexpressed
      ))

      df <- data.frame(
        target = t,
        sample = colnames(x),
        totalTags = x$totalTags,
        numberCTSSs = x$numberCTSSs,
        numberGenes = x$numberGenes
      )
      df[object$totalTags < t, c(
        "totalTags", "numberCTSSs",
        "numberGenes"
      )] <- NA

      df
    })
  )

  do.call("rbind", res)
}

#' Calculate the CTSS complexity of a CTSS dataset.
#'
#' This function calculates the CTSS complexity of a given CTSS
#' dataset.
#'
#' @param object A \code{RangedSummarizedExperiment} resulting from
#' CAGEfightR::quantifyCTSSs().
#' @param step The step size for the subsampling.
#' @param CTSSunexpressed The count threshold for considering a CTSS to be
#' unexpressed.
#' @param minCTSSsupport The minimum number of samples a CTSS must be
#' expressed in to be considered.
#'
#' @return A data frame containing the number of CTSSs (column "numberCTSSs")
#' for each sample and target (column "totalTags").
#'
#' @export
#'
#' @importFrom assertthat assert_that is.flag
#' @importFrom methods is
#'
calcCTSSComplexity <- function(object, step = 1e6, CTSSunexpressed = 1,
                               minCTSSsupport = 2) {
  assertthat::assert_that(
    methods::is(object, "RangedSummarizedExperiment"),
    is.numeric(step), length(step) == 1, step > 0,
    is.numeric(CTSSunexpressed), length(CTSSunexpressed) == 1, CTSSunexpressed >= 0,
    is.numeric(minCTSSsupport), length(minCTSSsupport) == 1, minCTSSsupport >= 1
  )
  object <- suppressWarnings(calcTotalTags(object,
    inputAssay = "counts"
  ))

  targets <- seq(step, max(object$totalTags), by = step)

  message("Subsampling counts")
  res <- c(
    list(data.frame(
      target = 0, sample = colnames(object),
      totalTags = 0, numberCTSSs = 0
    )),
    lapply(targets, function(t) {
      message(t)

      x <- subsampleTarget(object, "counts", t)
      x <- suppressWarnings(calcTotalTags(x, inputAssay = "counts"))
      if (minCTSSsupport > 1) {
        x <- suppressMessages(
          subsetBySupport(x,
            unexpressed = CTSSunexpressed,
            minSamples = minCTSSsupport
          )
        )
      }

      x <- suppressWarnings(calcNumberCTSSs(x,
        inputAssay = "counts",
        unexpressed = CTSSunexpressed
      ))

      df <- data.frame(
        target = t, sample = colnames(x),
        totalTags = x$totalTags, numberCTSSs = x$numberCTSSs
      )
      df[object$totalTags < t, c("totalTags", "numberCTSSs")] <- NA

      df
    })
  )

  do.call("rbind", res)
}

#' Calculate the gene complexity of a CTSS dataset.
#'
#' This function calculates the gene complexity of a given CTSS
#' dataset.
#'
#' @param object A \code{RangedSummarizedExperiment} resulting from
#' \code{CAGEfightR::quantifyCTSSs()}.
#' @param txModels A \code{TxDb} object containing the gene models.
#' @param step The step size for the subsampling.
#' @param CTSSunexpressed The count threshold for considering a CTSS to be
#' unexpressed.
#' @param geneunexpressed The count threshold for considering a gene to be
#' unexpressed.
#' @param minCTSSsupport The minimum number of samples a CTSS must be
#' expressed in to be considered.
#'
#' @return A data frame containing the number of genes (column "numberGenes")
#' for each sample and target (column "totalTags").
#'
#' @importFrom assertthat assert_that is.flag
#' @importFrom methods is
#'
#' @export
calcGeneComplexity <- function(object, txModels, step = 1e6,
                               CTSSunexpressed = 1, geneunexpressed = 9,
                               minCTSSsupport = 2) {
  assertthat::assert_that(
    methods::is(object, "RangedSummarizedExperiment"),
    methods::is(txModels, "TxDb"),
    is.numeric(step), length(step) == 1, step > 0,
    is.numeric(CTSSunexpressed), length(CTSSunexpressed) == 1, CTSSunexpressed >= 0,
    is.numeric(geneunexpressed), length(geneunexpressed) == 1, geneunexpressed >= 0,
    is.numeric(minCTSSsupport), length(minCTSSsupport) == 1, minCTSSsupport >= 1
  )
  object <- suppressWarnings(
    calcTotalTags(object, inputAssay = "counts")
  )
  object <- suppressMessages(
    suppressWarnings(
      assignGeneID(object,
        geneModels = txModels,
        outputColumn = "geneID"
      )
    )
  )

  targets <- seq(step, max(object$totalTags), by = step)

  message("Subsampling counts")
  res <- c(
    list(data.frame(
      target = 0, sample = colnames(object),
      totalTags = 0, numberGenes = 0
    )),
    lapply(targets, function(t) {
      message(t)

      x <- subsampleTarget(object, "counts", t)
      x <- suppressWarnings(
        calcTotalTags(x, inputAssay = "counts")
      )
      if (minCTSSsupport > 1) {
        x <- suppressMessages(
          subsetBySupport(x,
            unexpressed = CTSSunexpressed,
            minSamples = minCTSSsupport
          )
        )
      }

      x <- suppressWarnings(
        calcNumberGenes(x,
          txModels,
          inputAssay = "counts",
          unexpressed = geneunexpressed
        )
      )

      df <- data.frame(
        target = t, sample = colnames(x),
        totalTags = x$totalTags,
        numberGenes = x$numberGenes
      )
      df[object$totalTags < t, c("totalTags", "numberGenes")] <- NA

      df
    })
  )

  do.call("rbind", res)
}

#' Calculate the divergent loci complexity of a CTSS dataset.
#'
#' This function calculates the divergent loci complexity of a given CTSS
#' dataset.
#'
#' @param object A \code{RangedSummarizedExperiment} resulting from
#' CAGEfightR::quantifyCTSSs().
#' @param loci A \code{GRanges} object containing the tag clusters from which
#' divergent loci should be called.
#' @param step The step size for the subsampling.
#' @param CTSSunexpressed The count threshold for considering a CTSS to be
#' unexpressed.
#' @param lociunexpressed The count threshold for considering a tag cluster
#' to be unexpressed.
#' @param requirebidirectional Logical indicating whether to require
#' bidirectional transcription for a divergent locus from the same sample.
#' @param minCTSSsupport The minimum number of samples a CTSS must be
#' expressed in to be considered.
#'
#' @return A data frame containing the number of divergent loci (column
#' "numberLoci") for each sample and target (column "totalTags").
#' @importFrom assertthat assert_that is.flag
#' @importFrom methods is
#'
#' @export
calcDivergentLociComplexity <- function(object, loci, step = 1e6,
                                        CTSSunexpressed = 1,
                                        lociunexpressed = 2,
                                        requirebidirectional = FALSE,
                                        minCTSSsupport = 2) {
  assertthat::assert_that(
    methods::is(object, "RangedSummarizedExperiment"),
    methods::is(loci, "GRanges"),
    is.numeric(step), length(step) == 1, step > 0,
    is.numeric(CTSSunexpressed), length(CTSSunexpressed) == 1, CTSSunexpressed >= 0,
    is.numeric(lociunexpressed), length(lociunexpressed) == 1, lociunexpressed >= 0,
    assertthat::is.flag(requirebidirectional),
    is.numeric(minCTSSsupport), length(minCTSSsupport) == 1, minCTSSsupport >= 1
  )
  object <- suppressWarnings(
    calcTotalTags(object, inputAssay = "counts")
  )

  targets <- seq(step, max(object$totalTags), by = step)

  message("Subsampling counts")
  res <- c(
    list(data.frame(
      target = 0, sample = colnames(object),
      totalTags = 0, numberLoci = 0
    )),
    lapply(targets, function(t) {
      message(t)

      x <- subsampleTarget(object, "counts", t)
      x <- suppressWarnings(
        calcTotalTags(x, inputAssay = "counts")
      )
      if (minCTSSsupport > 1) {
        x <- suppressMessages(
          subsetBySupport(x,
            unexpressed = CTSSunexpressed,
            minSamples = minCTSSsupport
          )
        )
      }

      x <- suppressWarnings(
        calcNumberDivergentLoci(x,
          loci = loci, inputAssay = "counts",
          unexpressed = lociunexpressed,
          requirebidirectional = requirebidirectional
        )
      )

      df <- data.frame(
        target = t, sample = colnames(x),
        totalTags = x$totalTags, numberLoci = x$numberLoci
      )
      df[object$totalTags < t, c("totalTags", "numberLoci")] <- NA

      df
    })
  )

  do.call("rbind", res)
}

#' Calculate the complexity of a CTSS dataset and generate data for a fingerprint plot.
#'
#' This function calculates the complexity of a given CTSS dataset and generates a data
#' frame suitable for creating a fingerprint plot.
#' @param object A \code{RangedSummarizedExperiment} resulting from
#' \code{CAGEfightR::quantifyCTSSs()}.
#' @param replicates A character vector of replicate names to include in the analysis.
#' @param inputAssay A character string of the assay to use for the counts. Default is "counts".
#' @param n The number of points to include in the fingerprint plot. Default is 100,000.
#'
#' @return A data frame containing the sample names, raw counts, cumulative sums,
#' relative cumulative sums, and ranks for each sample.
#' @export
#'
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' @importFrom assertthat assert_that is.string is.numeric is.character
#' @importFrom methods is
#'
fingerPrint <- function(object, replicates, inputAssay = "counts", n = 100000) {
  assertthat::assert_that(
    methods::is(object, "RangedSummarizedExperiment"),
    assertthat::is.string(inputAssay),
    inputAssay %in% SummarizedExperiment::assayNames(object),
    is.character(replicates),
    length(replicates) >= 1,
    all(replicates %in% colnames(object)),
    is.numeric(n), length(n) == 1, n >= 1
  )

  all <- base::data.frame()
  for (i in 1:base::length(replicates)) {
    sam <- replicates[i]
    col <- base::sort(
      SummarizedExperiment::assay(object, inputAssay)[, sam]
    )

    ## Calculate cumSum as a function of CTSSs rank
    csum <- base::cumsum(col)
    totalLength <- base::length(csum)
    totalCTSSs <- csum[totalLength]
    csumRel <- csum / totalCTSSs
    rank <- base::seq(1L, totalLength, 1L) / totalLength

    ## Subset to memory-efficient 100k values
    subset <- base::c(
      base::seq(1, totalLength, base::floor(totalLength / n)),
      totalLength
    )

    df <- tibble::tibble(
      sample = sam,
      raw = col[subset],
      csum = csum[subset],
      csumRel = csumRel[subset],
      rank = rank[subset]
    )
    ## Combine
    all <- dplyr::bind_rows(all, df)
  }
  return(all)
}