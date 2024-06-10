#' Subsample a \code{SummarizedExperiment} to a target sequencing depth.
#' 
#' @param object A \code{SummarizedExperiment} object.
#' @param inputAssay The assay to use.
#' @param target The target sequencing depth.
#' 
#' @return A \code{SummarizedExperiment} object with the subsampled counts.
#' 
#' @export
#' 
#' @importFrom assertthat assert_that
#' @importFrom methods is
#' @import Matrix
#' @importFrom stats rbinom
#' 
subsampleTarget <- function(object, inputAssay = "counts", target) {
  
  assert_that(is(object, "SummarizedExperiment"),
                          inputAssay %in% assayNames(object),
                          is.numeric(target), target > 0)
  
  a <- assay(object,inputAssay)
  n <- ncol(a)
  nz <- lapply(1:n, function(i) nonzero(a[,i,drop=FALSE])[,1])
  s <- Matrix::colSums(a)
  d <- unlist(lapply(1:n, function(i) {
    if (s[i]<=target)
      a[nz[[i]],i]
    else
      rbinom(length(nz[[i]]),a[nz[[i]],i], target/s[i])
  }))
  keep <- which(sapply(nz,length)>0)
  assay(object, inputAssay) <- 
    Matrix::sparseMatrix(i=unlist(nz),
                         j=unlist(lapply(keep, function(i) {
                           rep(i,length(nz[[i]]))})),
                         x=d, dimnames=list(rownames(a), 
                                            colnames(a)[keep]))
  
  object
}

#' Subsample a \code{SummarizedExperiment} to a proportion of its sequencing 
#' depth. Useful for saturation analyses.
#' 
#' @param object A \code{SummarizedExperiment} object.
#' @param inputAssay The assay to use.
#' @param proportion The proportion of the sequencing depth to subsample.
#' 
#' @return A \code{SummarizedExperiment} object with the subsampled counts.
#' 
#' @export
subsampleProportion <- function(object, inputAssay = "counts", proportion) {
  
  assert_that(is(object, "SummarizedExperiment"),
                          inputAssay %in% assayNames(object),
                          is.numeric(proportion), 
                          proportion > 0 && proportion <=1)
  
  a <- assay(object,inputAssay)
  n <- ncol(a)
  nz <- lapply(1:n, function(i) nonzero(a[,i,drop=FALSE])[,1])
  d <- unlist(lapply(1:n, function(i) {
    rbinom(length(nz[[i]]),a[nz[[i]],i], proportion)
  }))
  keep <- which(sapply(nz,length)>0)
  assay(object, inputAssay) <- 
    Matrix::sparseMatrix(i=unlist(nz),
                         j=unlist(lapply(keep, function(i) {
                           rep(i,length(nz[[i]]))})),
                         x=d, dimnames=list(rownames(a), 
                                            colnames(a)[keep]))
  
  object
}
