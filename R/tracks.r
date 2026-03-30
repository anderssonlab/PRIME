#' Export bigWig files for each replicate, optionally merging strands.
#' 
#' @param object A \code{RangedSummarizedExperiment} object containing the CTSSs.
#' @param replicates A character vector of replicate names to export.
#' @param dir A character string of the directory to export the bigWig files to
#' @param inputAssay A character string of the assay to use for the bigWig scores. Default is "TMM".
#' @param splitByStrand A logical indicating whether to export strands separately. Default is FALSE.
#' 
#' @return No return value. The function exports bigWig files to the specified directory.
#' 
#' @importFrom assertthat assert_that
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom rtracklayer export
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment rowRanges
#' @export 
#' 
writeBw <- function(object, replicates, dir, inputAssay="TPM",splitByStrand = TRUE){
  
  ## Check consistency of replicates
  assertthat::assert_that(
      base::length(SummarizedExperiment::rowRanges(object)) == base::nrow(SummarizedExperiment::assay(object,inputAssay)),
      base::all(replicates %in% rownames(SummarizedExperiment::colData(object)))
    )
  base::stopifnot("Chosen directory doesn't exist!" = base::dir.exists(dir))
  
  gr <- GenomicRanges::GRanges(SummarizedExperiment::rowRanges(object))
  
  ## Parallel export across all chosen replicates
  foreach::foreach(i=1:base::length(replicates)) %dopar% {
    
    ## Print replicate when exiting
    base::on.exit(base::cat(base::paste0("Exit due to ",replicates[i],".\n")), add = TRUE)
    
    ## Convert individual replicate to GRanges object
    S4Vectors::mcols(gr)["score"] <- SummarizedExperiment::assay(object,inputAssay)[,replicates[i]]
    
    ## Export replicate's plus and minus strands separately
    if(splitByStrand){
      for (s in base::c("+","-")){
          grs <- gr[BiocGenerics::strand(gr) == s,]
          s <- base::ifelse(s == "+","plus","minus")
          
          fn <- base::paste0(dir,replicates[i],".",inputAssay,".",s,".bw")
          rtracklayer::export(
              object = grs,
              con = fn,
              format = "bigWig"
            )
          base::message(base::paste(fn,"done.",sep = " "))
        }
      } else {
        
        # Aggregate counts from both strands
        grred <- IRanges::reduce(gr, ignore.strand = TRUE, with.revmap = TRUE, min.gapwidth = 0)
        S4Vectors::mcols(grred) <- S4Vectors::aggregate(
            gr, S4Vectors::mcols(grred)$revmap, score = base::sum(score), drop = FALSE
          )
        
        # Keep only score metadata column
        grred@elementMetadata <- grred@elementMetadata[,"score",drop = FALSE]
    
        fn <- base::paste0(dir,replicates[i],".",inputAssay,".bw")
        rtracklayer::export(
              object = grred,
              con = fn,
              format = "bigWig"
            )
        base::message(base::paste(fn,"done.",sep = " "))
      }
    }
  }  