#' Calculate the distance to the median (DM) adjusted CV2 for a given 
#' \code{RangedSummarizedExperiment} object. Inspired by the method outlined in
#' Kolodziejczyk et al., 2015, Cell Stem Cell.
#' 
#' @param object A \code{RangedSummarizedExperiment} object.
#' @param inputAssay The name of the assay to use for the calculation.
#' @param prefix The prefix to use for the output columns.
#' 
#' @return A \code{RangedSummarizedExperiment} object with the adjusted CV2
#' in the metadata columns.
#' 
#' @export
#' 
#' @import SummarizedExperiment
#' @import Matrix
#' @importFrom zoo rollapply
#' 
DMadjustedCV <- function(object, inputAssay="TPM", prefix="") {
  
  object <- calcCV(object, inputAssay, paste(prefix,"CV",sep=""))
  object <- calcMean(object, inputAssay, paste(prefix,"mean",sep=""))
  mcols(object)[,paste(prefix,"log10_CV2",sep="")] <- 
    log10(rowData(object)[,paste(prefix,"CV",sep="")]^2)
  
  mean_order <- order(rowData(object)[,paste(prefix,"mean",sep="")])
  
  rolling_median <- rollapply(
    mcols(object)[,paste(prefix,"log10_CV2",sep="")][mean_order],
    width=50,by=25,FUN=median,
    fill=list("extend","extend","NA"))
  
  nas <- which(is.na(rolling_median))
  
  rolling_median[nas] <- median(
    mcols(object)[,paste(prefix,"log10_CV2",sep="")][mean_order][nas])
  
  names(rolling_median) <- (1:length(object))[mean_order]
  
  reorder <- match(1:length(object),names(rolling_median))
  
  rolling_median <- rolling_median[reorder]
  
  mcols(object)[, paste(prefix,"roll_median_log10_CV2",sep="")] <- 
    rolling_median
  
  mcols(object)[, paste(prefix,"adjusted_log10_CV2",sep="")] <- 
    mcols(object)[,paste(prefix,"log10_CV2",sep="")] - rolling_median
  
  object
}

## Helper function
calcCV <- function(object, inputAssay="counts", outputColumn="CV") {
  
  data <- assay(object, inputAssay)
  value <- rowSds(data,na.rm=TRUE)/Matrix::rowMeans(data,na.rm=TRUE)
  
  mcols(object)[, outputColumn] <- value
  
  object
}

## Helper function
calcMean <- function(object, inputAssay="counts", outputColumn="mean") {
  
  data <- assay(object, inputAssay)
  value <- Matrix::rowMeans(data,na.rm=TRUE)
  
  mcols(object)[, outputColumn] <- value
  
  object
}

## Helper function
rowVars <- function(x, ...) {
  sqr <- function(x)  x*x
  n <- Matrix::rowSums(!is.na(x))
  n[n<=1] <- NA
  return(rowSums(sqr(x-Matrix::rowMeans(x, ...)), ...)/(n-1))
}

## Helper function
rowSds <- function(x, ...)
  sqrt(rowVars(x, ...))


## Helper function
calcMedian <- function(object, inputAssay="counts", outputColumn="median") {
  
  data <- assay(object, inputAssay)
  value <- apply(data,1,median,na.rm=TRUE)
  
  mcols(object)[, outputColumn] <- value
  
  object
}

## Helper function
calcMax <- function(object, inputAssay="counts", outputColumn="median") {
  
  data <- assay(object, inputAssay)
  value <- apply(data,1,max,na.rm=TRUE)
  
  mcols(object)[, outputColumn] <- value
  
  object
}
