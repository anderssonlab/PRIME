
#' Extract heatmap-like CTSS data for regions of fixed width.
#' 
#' @param regions \code{GRanges} with regions to extract data for 
#' (must be of the same size).
#' @param data \code{GRanges} with value to consider in specified column 
#' (score).
#' @param column column in \code{mcols(data)} to consider.
#' @param transform_fn function applied on data for each region 
#' (e.g. to reduce dimension, defaults to identity).
#' @param ... additional arguments to pass to \code{transform_fn}.
#' 
#' @return list of matrices with data for each strand.
#' 
#' @export
#' 
#' @import GenomicRanges
#' @import IRanges
#' @import CAGEfightR
#' @importFrom assertthat assert_that
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' 
heatmapData <- function(regions, data, column="score", 
                        transform_fn=identity, ...) {
  
  ## Check regions
  assert_that(length(unique(width(regions)))==1,
                          column %in% colnames(mcols(data)))
  
  sl <- intersect(seqlevels(regions),seqlevels(data))
  if (!all(seqlevels(regions) %in% sl) || !all(seqlevels(data) %in% sl)) {
    warning(paste0("seqlevels differ between regions and data GRanges",
                   "objects, subsetting to intersection"))
    GenomeInfoDb::seqlevels(regions, pruning.mode="coarse") <- sl
    GenomeInfoDb::seqlevels(data, pruning.mode="coarse") <- sl
  }   
  
  regionsByStrand <- CAGEfightR:::splitByStrand(regions)
  dataByStrand <- splitPooledWeight(data, weight=column)
  
  nr <- names(regionsByStrand)[sapply(regionsByStrand, 
                                      function(x) length(x)>0)]
  nd <- names(dataByStrand)[sapply(dataByStrand, function(x) {
    any(sapply(x,function(y) any(y!=0)))})]
  
  res <- lapply(nr, function(r) {
    message("extracting data for strand: ", r)
    
    dat <- lapply(nd, function(d) {
      message("   ",d)
      vl <- Views(dataByStrand[[d]], regionsByStrand[[r]])
      n <- paste0(unlist(lapply(names(vl), function(nv) {
        rep(nv,length(vl[[nv]]))})), ":",
        unlist(lapply(vl,start)), "-", 
        unlist(lapply(vl,end)), ";", d)
      
      vlen <- sapply(vl,length)
      
      m <- do.call("rbind",
                   lapply(vl[vlen>0], function(v) {
                     vectorListToMatrix( 
                       viewApply(v,function(x) {
                         as(transform_fn(as.vector(x,mode="numeric"),...),"dsparseVector")
                       })
                     )}))
      
      rownames(m) <- n
      m
    })
    names(dat) <- nd
    
    dat
  })
  
  names(res) <- nr
  
  res
}

## Helper function, not exported. Converts a list of sparse vectors to a sparse matrix.
vectorListToMatrix <- function(vl){  
  sm_i<-NULL
  sm_j<-NULL
  sm_x<-NULL
  for (k in 1:length(vl)) {
    sm_i <- c(sm_i,rep(k,length(vl[[k]]@i)))
    sm_j <- c(sm_j,vl[[k]]@i)
    sm_x <- c(sm_x,vl[[k]]@x)
  }
  return (sparseMatrix(i=sm_i,j=sm_j,x=sm_x,dims=c(length(vl),vl[[1]]@length)))
}

## Helper function, not exported, adds functionality to CAGEfightR:::splitPooled
## by allowing to focus on another column than "score"
splitPooledWeight <- function(object, weight="score"){

    ## Split by strand
    o <- CAGEfightR:::splitByStrand(object)

    ## Calculate coverage
    o <- lapply(o, coverage, weight=weight)

    ## Round to handle floating point errors
    o <- lapply(o, round, digits=9)

    ## Return
    o
}
        
