
#' Conditional normalization method (e.g. on GC content) based on approach 
#' described in Pickrell, et al, 2010, Nature.
#' 
#' @param object A \code{SummarizedExperiment} object.
#' @param inputAssay The name of the assay to use for the calculation.
#' @param outputAssay The name of the assay to use for the output.
#' @param conditionalColumn The column in the metadata to use for the 
#' conditional normalization.
#' @param offsetAssay The name of the assay to use for the offset.
#' @param bins The number of bins to use for the conditional normalization.
#' @param sizeFactors The size factors to use for the normalization.
#' @param minCount The minimum count to consider a unit (CTSS, gene, etc.) 
#' as expressed.
#' @param minSamples The minimum number of samples to consider a unit 
#' as expressed.
#' @param aggregate.fn The function to use for aggregation.
#' 
#' @return A \code{SummarizedExperiment} object with the normalized counts.
#' 
#' @export
#' 
#' @importFrom Hmisc cut2
#' @import S4Vectors
#' @importFrom stats smooth.spline predict
#' @importFrom assertthat assert_that
#' @importFrom methods is
#' 
conditionalNormalize <- function(object, inputAssay="counts", 
                                 outputAssay="normalized", 
                                 conditionalColumn="GC", offsetAssay=NULL, 
                                 bins=200, sizeFactors=NULL, minCount=1, 
                                 minSamples=1,aggregate.fn=sum) {
  
  assert_that(is(object, "SummarizedExperiment"),
                          inputAssay %in% assayNames(object),
                          conditionalColumn %in% colnames(mcols(object)))
  
  message("binning data according to conditional column...")
  
  y <- assay(object,inputAssay)
  keep <- which(Matrix::rowSums(y>minCount) >= minSamples)
  
  x <- as.numeric(mcols(object)[,conditionalColumn])
  b <- as.character(Hmisc::cut2(x, unique(quantile(x, seq(1/bins, 1, 1/bins)))))
  b.m <- aggregate(x, by=list(b), FUN=mean)
  b.m.keep <- aggregate(x[keep], by=list(b[keep]), FUN=mean)
  n <- b.m[,1]
  b.m <- b.m[,2]
  names(b.m) <- n
  n <- b.m.keep[,1]
  b.m.keep <- b.m.keep[,2]
  names(b.m.keep) <- n
  
  y.b <- aggregate(as.matrix(y[keep,]),by=list(b[keep]),FUN=aggregate.fn)
  b.g <- y.b[,1]
  rownames(y.b) <- b.g
  y.b <- y.b[names(b.m.keep),]
  y.b <- data.matrix(y.b[,-1])
  
  message("calculating relative enrichments...")
  
  if (is.null(sizeFactors))
    sizeFactors <- rep(1,ncol(y))
  else
    sizeFactors <- sizeFactors / mean(sizeFactors)
  
  y.b <- sapply(1:ncol(y.b), function(i) y.b[,i] / sizeFactors[i])
  s <- colSums(y.b) / sum(y.b)
  
  f <- log2( t( t(y.b / rowSums(y.b)) / s ) )
  fit <- lapply(1:ncol(y), function(i) {
    b.m.i <- b.m.keep
    f.i <- f[,i]
    missing <- which(is.na(f[,i]) | is.infinite(f[,i]))
    if (length(missing)>0) {
      b.m.i <- b.m.i[-missing]
      f.i <- f.i[-missing]
    }
    stats::smooth.spline(b.m.i, f.i, spar=1)
  })
  
  message("normalizing data according to enrichments...")
  
  offset <- -1*Matrix::Matrix(sapply(1:ncol(y), function(i) {
    stats::predict(fit[[i]],b.m[b])$y}))
  y_normed <- Matrix::Matrix(round(y * 2^(offset)))
  y_normed <- Matrix::Matrix(sapply(1:ncol(y), function(i) {
    x <- y_normed[,i]
    x[x==0 & y[,i]>0] <- 1; x}))
  
  dimnames(y_normed) <- dimnames(y)
  dimnames(offset) <- dimnames(y)
  
  assay(object, outputAssay) <- y_normed
  if (!is.null(offsetAssay))
    assay(object, offsetAssay) <- offset
  
  ## return
  object
}

#' Normalize expression counts by size factors.
#' 
#' @param object A \code{SummarizedExperiment} object.
#' @param sizeFactors The size factors to use for the normalization.
#' @param inputAssay The assay to use.
#' @param outputAssay The assay to use for the output.
#' 
#' @return A \code{SummarizedExperiment} object with the normalized counts.
#' 
#' @export
#' 
#' @importFrom Matrix t rowSums colSums
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom assertthat assert_that
#' @importFrom methods is
#' 
normalizeBySizeFactors <- function(object, sizeFactors, inputAssay="counts", 
                                   outputAssay="normalized") {
  
  assert_that(is(object, "SummarizedExperiment"),
                          inputAssay %in% assayNames(object),
                          length(sizeFactors) == ncol(object))
  
  assay(object, outputAssay) <- Matrix::t(Matrix::t(
    assay(object, inputAssay)) / sizeFactors)
  
  object
}

#' TPM normalize expression counts, considering size factors.
#' 
#' @param object A \code{SummarizedExperiment} object.
#' @param sizeFactors The size factors to use for the normalization.
#' @param inputAssay The assay to use.
#' @param outputAssay The assay to use for the output.
#' 
#' @return A \code{SummarizedExperiment} object with the TPM normalized counts.
#' 
#' @export
TPMnormalizeBySizeFactors <- function(object, sizeFactors, inputAssay="counts", 
                                      outputAssay="TPM") {
  
  assert_that(is(object, "SummarizedExperiment"),
                          inputAssay %in% assayNames(object),
                          length(sizeFactors) == ncol(object))
  
  ## Centre size factors
  sf <- sizeFactors / mean(sizeFactors)
  
  ## Calculate scaling factors
  a <- assay(object, inputAssay)
  sf <- 1e6 / (sf * mean(Matrix::colSums(a)))
  
  assay(object, outputAssay) <- Matrix::t(Matrix::t(a) * sf)
  
  object
}
