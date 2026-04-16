
## Helper function
## Calculate number of unique CTSS positions
#' @importFrom assertthat assert_that not_empty is.string
#' @importFrom methods is
#' @importFrom CAGEfightR calcTotalTags assignGeneID quantifyGenes
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom S4Vectors mcols
#' @importFrom Matrix colSums

calcNumberCTSSs <- function(object, inputAssay = "counts", 
                            outputColumn = "numberCTSSs", 
                            unexpressed = 0) {
  ## Prechecks
  assertthat::assert_that(
    methods::is(object, "SummarizedExperiment"),
    assertthat::not_empty(object),
    inputAssay %in% SummarizedExperiment::assayNames(object),
    assertthat::is.string(inputAssay),
    assertthat::is.string(outputColumn)
  )
  
  if (outputColumn %in% colnames(SummarizedExperiment::colData(object))) {
    warning("object already has a column named ", outputColumn, 
            " in colData: It will be overwritten!")
  }
  
  ## Calculate number of unique CTSS positions
  a <- SummarizedExperiment::assay(object, inputAssay)
  SummarizedExperiment::colData(object)[, outputColumn] <- sapply(colnames(a), 
                                            function(i) sum(a[,i]>unexpressed))
  
  ## Return
  object
}

## Helper function
## Calculate number of unique genes
calcNumberGenes <- function(object, txModels, inputAssay = "counts", 
                            outputColumn = "numberGenes", 
                            unexpressed = 0) {
  ## Prechecks
  assertthat::assert_that(
    methods::is(object, "SummarizedExperiment"),
    assertthat::not_empty(object),
    inputAssay %in% SummarizedExperiment::assayNames(object),
    assertthat::is.string(inputAssay),
    assertthat::is.string(outputColumn)
  )
  
  if (outputColumn %in% colnames(SummarizedExperiment::colData(object))) {
    warning("object already has a column named ", outputColumn, 
            " in colData: It will be overwritten!")
  }
  
  if (!"geneID" %in% colnames(S4Vectors::mcols(object)))
    object <- CAGEfightR::assignGeneID(object, geneModels = txModels, 
                           outputColumn = "geneID")
  genelevel <- CAGEfightR::quantifyGenes(object, genes="geneID", 
                                         inputAssay=inputAssay)
  
  ## Calculate number of expressed genes
  a <- SummarizedExperiment::assay(genelevel, inputAssay)
  SummarizedExperiment::colData(object)[, outputColumn] <- sapply(colnames(a), 
                                            function(i) sum(a[,i]>unexpressed))
  
  ## Return
  object
}

## Helper function
## Calculate number of divergent loci
calcNumberDivergentLoci <- function(object, loci, inputAssay="counts", 
                                    outputColumn = "numberLoci", 
                                    unexpressed = 0, 
                                    requirebidirectional=FALSE) {
  ## Prechecks
  assertthat::assert_that(
    methods::is(object, "SummarizedExperiment"),
    assertthat::not_empty(object),
    inputAssay %in% SummarizedExperiment::assayNames(object),
    assertthat::is.string(inputAssay),
    assertthat::is.string(outputColumn)
  )
  
  if (outputColumn %in% colnames(SummarizedExperiment::colData(object))) {
    warning("object already has a column named ", outputColumn, 
            " in colData: It will be overwritten!")
  }
  
  ## Calculate number of expressed loci
  res <- quantifyStrandwiseDivergentLoci(loci, object, inputAssay=inputAssay)
  m <- SummarizedExperiment::assay(res$'-', inputAssay)
  p <- SummarizedExperiment::assay(res$'+', inputAssay)
  expressed <- sapply(colnames(m), function(i) m[,i]+p[,i]>unexpressed)
  if (requirebidirectional)
    expressed <- expressed & sapply(colnames(m), 
                                    function(i) m[,i]>0 & p[,i]>0)
  SummarizedExperiment::colData(object)[, outputColumn] <- Matrix::colSums(expressed)
  
  ## Return
  object
}