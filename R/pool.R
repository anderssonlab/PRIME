#' Pools expression signal across replicates
#' 
#' @param object  A \code{SummarizedExperiment} object.
#' @param replicates A \code{factor} indicating the replicate ID of each sample.
#' @param inputAssay The assay to use
#' 
#' @return A \code{SummarizedExperiment} expression pooled across replicates and merged column information
#' 
#' @export 
#' 
#' @importFrom Matrix rowSums
#' @importFrom assertthat assert_that
#' @import SummarizedExperiment
#' 

poolReplicates <- function(object, replicates, inputAssay="counts") {
    
    ## dimensions
    assert_that(ncol(object)==length(replicates),
                length(unique(replicates))>1)


    message("Pooling counts...")
    a <- do.call("cbind", lapply(unique(replicates), function(b) 
      Matrix::rowSums(assay(object,inputAssay)[,replicates==b,drop=FALSE])))
    a <- as(a, "sparseMatrix")
    colnames(a) <- unique(replicates)

    message("Create a pasted design matrix...")
    c <- do.call("rbind", lapply(unique(replicates), function(b) 
      apply(colData(object)[replicates==b,], 2 , paste , collapse = "," )))

    rse <- SummarizedExperiment(assays=SimpleList(counts=a),
                            rowRanges=rowRanges(object), colData=as.data.frame(c))
    
    rse <- calcTotalTags(rse)
    
    rse
  }
