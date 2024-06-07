#' Calculate the gene and CTSS complexity of a CTSS dataset.
#'
#' This function calculates the CTSS and gene complexity of a given CTSS
#' dataset.
#'
#' @param object A RangedSummarizedExperiment resulting from
#' CAGEfightR::quantifyCTSSs().
#' @param txModels A TxDb object containing the gene models.
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
calcComplexity <- function(object, txModels, step = 1e6, CTSSunexpressed = 1,
                           geneunexpressed = 9, minCTSSsupport = 2) {
  
  object <- suppressWarnings(CAGEfightR::calcTotalTags(object, 
                                                       inputAssay = "counts"))
  object <- suppressMessages(suppressWarnings(
    CAGEfightR::assignGeneID(object, 
                             geneModels = txModels, 
                             outputColumn = "geneID")))
  
  targets <- seq(step, max(object$totalTags), by = step)
  
  message("Subsampling counts")
  res <- c(
    list(data.frame(target = 0, sample = colnames(object), totalTags = 0, 
                    numberCTSSs = 0, numberGenes = 0)),
    lapply(targets, function(t) {
      message(t)
      
      x <- subsampleTarget(object, "counts", t)
      x <- suppressWarnings(CAGEfightR::calcTotalTags(x, inputAssay = "counts"))
      if (minCTSSsupport > 1) {
        x <- suppressMessages(CAGEfightR::subsetBySupport(x, 
                                              unexpressed = CTSSunexpressed, 
                                              minSamples = minCTSSsupport))
      }
      
      x <- suppressWarnings(calcNumberCTSSs(x, 
                                            inputAssay = "counts", 
                                            unexpressed = CTSSunexpressed))
      x <- suppressWarnings(calcNumberGenes(x, 
                                            txModels, 
                                            inputAssay = "counts", 
                                            unexpressed = geneunexpressed))
      
      df <- data.frame(target = t, 
                       sample = colnames(x), 
                       totalTags = x$totalTags, 
                       numberCTSSs = x$numberCTSSs, 
                       numberGenes = x$numberGenes)
      df[object$totalTags < t, c("totalTags", "numberCTSSs", 
                                 "numberGenes")] <- NA
      
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
#' @param object A RangedSummarizedExperiment resulting from
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
calcCTSSComplexity <- function(object, step = 1e6, CTSSunexpressed = 1, 
                               minCTSSsupport = 2) {
  
  object <- suppressWarnings(CAGEfightR::calcTotalTags(object, 
                                                       inputAssay = "counts"))
  
  targets <- seq(step, max(object$totalTags), by = step)
  
  message("Subsampling counts")
  res <- c(
    list(data.frame(target = 0, sample = colnames(object), 
                    totalTags = 0, numberCTSSs = 0)),
    lapply(targets, function(t) {
      message(t)
      
      x <- subsampleTarget(object, "counts", t)
      x <- suppressWarnings(CAGEfightR::calcTotalTags(x, inputAssay = "counts"))
      if (minCTSSsupport > 1) {
        x <- suppressMessages(
          CAGEfightR::subsetBySupport(x, 
                                      unexpressed = CTSSunexpressed, 
                                      minSamples = minCTSSsupport))
      }
      
      x <- suppressWarnings(calcNumberCTSSs(x, 
                                            inputAssay = "counts", 
                                            unexpressed = CTSSunexpressed))
      
      df <- data.frame(target = t, sample = colnames(x), 
                       totalTags = x$totalTags, numberCTSSs = x$numberCTSSs)
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
#' @param object A RangedSummarizedExperiment resulting from
#' CAGEfightR::quantifyCTSSs().
#' @param txModels A TxDb object containing the gene models.
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
#' @export
calcGeneComplexity <- function(object, txModels, step = 1e6, 
                               CTSSunexpressed = 1, geneunexpressed = 9, 
                               minCTSSsupport = 2) {
  object <- suppressWarnings(
    CAGEfightR::calcTotalTags(object, inputAssay = "counts"))
  object <- suppressMessages(
    suppressWarnings(
      CAGEfightR::assignGeneID(object, 
                               geneModels = txModels, 
                               outputColumn = "geneID")))
  
  targets <- seq(step, max(object$totalTags), by = step)
  
  message("Subsampling counts")
  res <- c(
    list(data.frame(target = 0, sample = colnames(object), 
                    totalTags = 0, numberGenes = 0)),
    lapply(targets, function(t) {
      message(t)
      
      x <- subsampleTarget(object, "counts", t)
      x <- suppressWarnings(
        CAGEfightR::calcTotalTags(x, inputAssay = "counts"))
      if (minCTSSsupport > 1) {
        x <- suppressMessages(
          CAGEfightR::subsetBySupport(x, 
                                      unexpressed = CTSSunexpressed, 
                                      minSamples = minCTSSsupport))
      }
      
      x <- suppressWarnings(
        CAGEfightR::calcNumberGenes(x, 
                                    txModels, 
                                    inputAssay = "counts", 
                                    unexpressed = geneunexpressed))
      
      df <- data.frame(target = t, sample = colnames(x), 
                       totalTags = x$totalTags, 
                       numberGenes = x$numberGenes)
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
#' @param object A RangedSummarizedExperiment resulting from
#' CAGEfightR::quantifyCTSSs().
#' @param loci A GRanges object containing the tag clusters from which 
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
#'
#' @export
calcDivergentLociComplexity <- function(object, loci, step = 1e6, 
                                        CTSSunexpressed = 1, 
                                        lociunexpressed = 2, 
                                        requirebidirectional = FALSE, 
                                        minCTSSsupport = 2) {
  
  object <- suppressWarnings(
    CAGEfightR::calcTotalTags(object, inputAssay = "counts"))
  
  targets <- seq(step, max(object$totalTags), by = step)
  
  message("Subsampling counts")
  res <- c(
    list(data.frame(target = 0, sample = colnames(object), 
                    totalTags = 0, numberLoci = 0)),
    lapply(targets, function(t) {
      message(t)
      
      x <- subsampleTarget(object, "counts", t)
      x <- suppressWarnings(
        CAGEfightR::calcTotalTags(x, inputAssay = "counts"))
      if (minCTSSsupport > 1) {
        x <- suppressMessages(
          CAGEfightR::subsetBySupport(x, unexpressed = CTSSunexpressed, 
                                      minSamples = minCTSSsupport))
      }
      
      x <- suppressWarnings(
        calcNumberDivergentLoci(x, loci = loci, inputAssay = "counts", 
                                unexpressed = lociunexpressed, 
                                requirebidirectional = requirebidirectional))
      
      df <- data.frame(target = t, sample = colnames(x), 
                       totalTags = x$totalTags, numberLoci = x$numberLoci)
      df[object$totalTags < t, c("totalTags", "numberLoci")] <- NA
      
      df
    })
  )
  
  do.call("rbind", res)
}
