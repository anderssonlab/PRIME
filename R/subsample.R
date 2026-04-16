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
#' @importFrom Matrix colSums sparseMatrix
#' @importFrom stats rbinom
#' 
subsampleTarget <- function(object, inputAssay = "counts", target) {

    assertthat::assert_that(
        methods::is(object, "SummarizedExperiment"),
        inputAssay %in% SummarizedExperiment::assayNames(object),
        base::is.numeric(target),
        target > 0
      )
    
    a <- SummarizedExperiment::assay(object,inputAssay)
    n <- base::ncol(a)
    nz <- base::lapply(
      X = 1:n,
      FUN = function(i){
          nonzero(a[,i,drop=FALSE])[,1]
        }
      )
    s <- Matrix::colSums(a)
    d <- base::unlist(
      base::lapply(
        X = 1:n,
        FUN = function(i){
        if (s[i]<=target){
              a[nz[[i]],i]
            } else {
              stats::rbinom(base::length(nz[[i]]),a[nz[[i]],i], target/s[i])
            }
          }
        )
      )
    keep <- base::which(
        base::sapply(nz,length)>0
      )
    
    # Create a matrix to store the subsampled counts, ensuring all rows are retained
    new_matrix <- Matrix::sparseMatrix(
      i = base::unlist(nz),
      j = base::unlist(
        base::lapply(
          X = keep,
          FUN = function(i){
              base::rep(i, base::length(nz[[i]]))
            }
          )
        ),
      x = d, 
      dims = base::dim(a),  # Preserve the original dimensions
      dimnames = base::list(
          base::rownames(a),
          base::colnames(a)[keep]
        )
      )
    
    # Replace the original assay with the new matrix
    SummarizedExperiment::assay(object, inputAssay) <- new_matrix
    
    # Recalculate the total count of reads for each sample
    object <- CAGEfightR::calcTotalTags(object)
    
    return(object)
  }


#' Subsample a \code{SummarizedExperiment} to a proportion of its sequencing 
#' depth. Useful for saturation analyses.
#' 
#' @param object A \code{SummarizedExperiment} object.
#' @param inputAssay The assay to use.
#' @param proportion The proportion of the sequencing depth to subsample.
#' 
#' @importFrom assertthat assert_that
#' @importFrom methods is
#' @importFrom Matrix sparseMatrix
#' @importFrom CAGEfightR calcTotalTags
#' @importFrom stats rbinom
#' @importFrom SummarizedExperiment assay assayNames
#'
#' @export
subsampleProportion <- function(object, inputAssay = "counts", proportion) {
  
  assertthat::assert_that(
    methods::is(object, "SummarizedExperiment"),
    inputAssay %in% SummarizedExperiment::assayNames(object),
    is.numeric(proportion), 
    proportion > 0 && proportion <=1
  )
  
  a <- SummarizedExperiment::assay(object,inputAssay)
  n <- ncol(a)
  ##nz <- lapply(1:n, function(i) DAPAR::nonzero(a[,i,drop=FALSE])[,1])
  nz <- lapply(1:n, function(i) nonzero(a[,i,drop=FALSE])[,1])
  d <- unlist(lapply(1:n, function(i) {
    stats::rbinom(length(nz[[i]]),a[nz[[i]],i], proportion)
  }))
  keep <- which(sapply(nz,length)>0)
  SummarizedExperiment::assay(object, inputAssay) <- 
    Matrix::sparseMatrix(i=unlist(nz),
                         j=unlist(lapply(keep, function(i) {
                           rep(i,length(nz[[i]]))})),
                         x=d, 
                         dims = dim(a),  # Preserve the original dimensions
                         dimnames=list(rownames(a), 
                         colnames(a)[keep]))

  # Recalculate the total count of reads for each sample
  object <- CAGEfightR::calcTotalTags(object)
                      
  return(object)
}

#' nonzero function from DAPAR package (https://github.com/edyp-lab/DAPAR)
#' temporarily including a local version to avoid installation issues
#' 
#' @param x A sparse matrix of class "dgCMatrix".
#' @return A two-column matrix containing the indices of the non-zero elements in the input matrix.
#' 
nonzero <- function(x) {
    ## function to get a two-column matrix containing the indices of the
    ### non-zero elements in a "dgCMatrix" class matrix

    stopifnot(inherits(x, "dgCMatrix"))
    if (all(x@p == 0)) {
        return(matrix(0,
            nrow = 0, ncol = 2,
            dimnames = list(character(0), c("row", "col"))
        ))
    }
    res <- cbind(x@i + 1, rep(seq(dim(x)[2]), diff(x@p)))
    colnames(res) <- c("row", "col")
    res <- res[x@x != 0, , drop = FALSE]
    return(res)
}
