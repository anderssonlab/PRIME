#' Initiator classification of \code{GRanges} into YR, YC and other.
#' 
#' @param object A \code{GRanges} or \code{StitchedGPos} object.
#' @param bsg The \code{BSgenome} reference to use. 
#' 
#' @return A \code{GRanges} object initiator classification.
#' 
#' @export
#' 
#' @import SummarizedExperiment
#' @import IRanges
#' @import CAGEfightR
#' @import BSgenome
#' 
#' 

initiatorScore <- function(object, bsg) {
  # Ensure the input is of class GRanges or StitchedGPos
  if (!inherits(object, c("GRanges", "StitchedGPos"))) {
    stop("object must be of class GRanges or StitchedGPos")
  }
  
  # Remove out-of-bound indices
  idx <- GenomicRanges:::get_out_of_bound_index(object)
  if (length(idx) > 0) {
    object <- object[-idx]
  }
  
  # Create promoters based on the class of object
  if (class(object) == "GRanges") {
    swapped_object <- swapRanges(object)
    x <- promoters(swapped_object, upstream = 1, downstream = 1)
  } else if (class(object)== "StitchedGPos") {
    x <- promoters(object, upstream = 1, downstream = 1)
  }
    
  # Get sequences
  x$seq <- getSeq(bsg, x, as.character = TRUE)
  
  # Define INR categories
  YC <- c("CC", "TC")
  YR <- c("CA", "TA", "CG", "TG")
  
  # Assign INR category
  x$INR <- "other"
  x$INR[x$seq %in% YC] <- "YC"
  x$INR[x$seq %in% YR] <- "YR"
  
  return(x)
}