#' Export bigWig files for each replicate.
#' 
#' @param object A \code{RangedSummarizedExperiment} object containing the CTSSs.
#' @param replicates A character vector of replicate names to export.
#' @param dir A character string of the directory to export the bigWig files to.
#' @param inputAssay A character string of the assay to use for the bigWig scores. Default is "TMM".
#' 
#' @return No return value. The function exports bigWig files to the specified directory.
#' 
#' @import assertthat
#' @import foreach
#' @import rtracklayer
#' @import SummarizedExperiment
#' 
#' @export
#' 

writeBw <- function(object, replicates, dir, inputAssay="TMM"){
  
  ## Check consistency of replicates
  assertthat::assert_that(
      base::length(SummarizedExperiment::rowRanges(object)) == base::nrow(SummarizedExperiment::assay(object,"TMM")),
      base::all(replicates %in% rownames(SummarizedExperiment::colData(object)))
    )
  base::stopifnot("Chosen directory doesn't exist!" = base::dir.exists(dir))

  i <- NULL  # avoid R CMD check NOTE
  
  ## Parallel export across all chosen replicates
  foreach::foreach(i=1:base::length(replicates)) %dopar% {
    
    ## Convert individual replicate to GRanges object
    gr <- SummarizedExperiment::rowRanges(object)
    S4Vectors::mcols(gr)["score"] <- SummarizedExperiment::assay(object,inputAssay)[,replicates[i]]
    
    ## Export replicate's plus and minus strands separately
    for (s in base::c("+","-")){
        grs <- gr[BiocGenerics::strand(gr) == s,]
        s <- base::ifelse(s == "+","plus","minus")
        fn <- base::paste0(dir,replicates[i],".tmm.",s,".bw")
        
        rtracklayer::export(
            object = grs,
            con = fn,
            format = "bigWig"
          )
        base::message(base::paste(fn,"done.",sep = " "))
      }
    }
  }