#' Calculate batch support for expression data.
#' 
#' @param object A \code{SummarizedExperiment} object.
#' @param batch A \code{factor} indicating the batch of each sample.
#' @param inputAssay The assay to use.
#' @param outputColumn The name of the column to store the batch support.
#' @param unexpressed The threshold below which a unit (CTSS, gene, etc.) 
#' is considered unexpressed.
#' @param minSamples The minimum number of samples in which a unit must be
#' expressed to be considered supported.
#' 
#' @return A \code{SummarizedExperiment} object with the batch support column.
#' 
#' @export
#' 
#' @importFrom Matrix rowSums
#' @importFrom assertthat assert_that
#' @import SummarizedExperiment
#' 
calcBatchSupport <- function(object, batch, inputAssay="counts", 
                             outputColumn="batchSupport", unexpressed=0, 
                             minSamples=1) {
  
  ## dimensions
  assert_that(ncol(object)==length(batch),
              length(unique(batch))>1,
              minSamples>=1)
  
  if (outputColumn %in% colnames(rowData(object))) {
    warning("object already has a column named ", outputColumn, 
            " in rowData: It will be overwritten!")
  }
  
  message("Calculating support per batch...")
  a <- 
    do.call("cbind", lapply(unique(batch), function(b) {
      Matrix::rowSums(assay(object,inputAssay)[,batch==b,drop=FALSE]>
                        unexpressed)
    }))
  
  message("Calculating overall batch support...")
  rowData(object)[, outputColumn] <- as.integer(Matrix::rowSums(a>=minSamples))
  
  object
}

#' Remove singleton counts from a \code{SummarizedExperiment} object. 
#' Singleton counts are those that have a value of 1 in the specified assay. 
#' This function replaces those values with 0 and stores the modified assay in a new column.
#' 
#' @param rse A \code{SummarizedExperiment} object.
#' @param inputAssay The name of the assay to use.
#' @param outputAssay The name of the assay to store the modified counts. Default is "counts.noSingletons".
#' 
#' @return A \code{SummarizedExperiment} object with the modified assay.
#' @export
#' 
rmSingletons <- function(rse, inputAssay = "counts", outputAssay = "counts.noSingletons") {
  rmS <- assay(rse, inputAssay)
  # Find the indices where values are equal to 1
  idx <- which(rmS == 1, arr.ind = TRUE)
  # Replace the values at those indices with 0
  rmS[idx] <- 0
  # Assign the modified assay back to the SummarizedExperiment object
  assay(rse, outputAssay) <- rmS
  rse
}
