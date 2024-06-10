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
