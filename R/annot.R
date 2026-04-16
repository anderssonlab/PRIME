#' Calculate the total tags for each annotation type.
#'
#' @param data A \code{RangedSummarizedExperiment} object containing the CTSSs.
#' @param txModels A \code{TxDb} object containing the the transcript models for annotation.
#' @param uniqueCTSS A \code{logical} indicating whether to count only unique CTSSs (TRUE) or to sum the total tags (FALSE).
#'
#' @return A data frame containing the total tags for each annotation type.
#' @export
#'
#' @importFrom CAGEfightR assignTxType calcTotalTags
#' @importFrom SummarizedExperiment assay assayNames colData rowData
#' @importFrom Matrix colSums
#' @importFrom assertthat assert_that is.flag
#' @importFrom methods is
#'

calcAnnoCTSS <- function(data, txModels, uniqueCTSS = FALSE) {
  assertthat::assert_that(
    methods::is(data, "RangedSummarizedExperiment"),
    methods::is(txModels, "TxDb"),
    assertthat::is.flag(uniqueCTSS)
  )

  annos <- CAGEfightR::assignTxType(
    object = data,
    txModels = txModels
  )

  assertthat::assert_that("txType" %in% colnames(SummarizedExperiment::rowData(annos)))

  # Identify all txTypes
  types <- base::unique(SummarizedExperiment::rowData(annos)$txType)

  # Initialize list to store total tags for each txType
  totalTags_list <- base::vector("list", length = base::length(types))

  # Loop through each txType
  for (i in base::seq_along(types)) {
    # Subset data based on current txType
    idx <- SummarizedExperiment::rowData(annos)$txType %in% types[i]
    current_subset <- annos[idx, ]
    # current_subset <- S4Vectors::subset(annos, txType %in% types[i])
    if (uniqueCTSS) {
      # Calculate totalTags for current subset
      current_subset <- Matrix::colSums(
        SummarizedExperiment::assay(current_subset) > 0
      )
      totalTags_list[[i]] <- current_subset
    } else {
      # Calculate totalTags for current subset
      current_subset <- SummarizedExperiment::colData(
        CAGEfightR::calcTotalTags(current_subset)
      )
      totalTags_list[[i]] <- current_subset$totalTags
    }
  }

  # Combine totalTags into a data frame
  ann_counts <- base::as.data.frame(
    base::do.call(base::cbind, totalTags_list)
  )
  base::colnames(ann_counts) <- types
  return(ann_counts)
}