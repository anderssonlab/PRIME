
## Helper function
## Calculate number of unique CTSS positions
#' @importFrom assertthat assert_that not_empty is.string
#' @importFrom methods is
#' @import CAGEfightR
#' @import SummarizedExperiment

calcNumberCTSSs <- function(object, inputAssay = "counts", 
                            outputColumn = "numberCTSSs", 
                            unexpressed = 0) {
  ## Prechecks
  assert_that(is(object, "SummarizedExperiment"),
                          not_empty(object),
                          inputAssay %in% assayNames(object),
                          is.string(inputAssay),
                          is.string(outputColumn))
  
  if (outputColumn %in% colnames(colData(object))) {
    warning("object already has a column named ", outputColumn, 
            " in colData: It will be overwritten!")
  }
  
  ## Calculate number of unique CTSS positions
  a <- assay(object, inputAssay)
  colData(object)[, outputColumn] <- sapply(colnames(a), 
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
  assert_that(is(object, "SummarizedExperiment"),
                          not_empty(object),
                          inputAssay %in% assayNames(object),
                          is.string(inputAssay),
                          is.string(outputColumn))
  
  if (outputColumn %in% colnames(colData(object))) {
    warning("object already has a column named ", outputColumn, 
            " in colData: It will be overwritten!")
  }
  
  if (!"geneID" %in% colnames(mcols(object)))
    object <- assignGeneID(object, geneModels = txModels, 
                           outputColumn = "geneID")
  genelevel <- quantifyGenes(object, genes="geneID", 
                                         inputAssay=inputAssay)
  
  ## Calculate number of expressed genes
  a <- assay(genelevel, inputAssay)
  colData(object)[, outputColumn] <- sapply(colnames(a), 
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
  assert_that(is(object, "SummarizedExperiment"),
                          not_empty(object),
                          inputAssay %in% assayNames(object),
                          is.string(inputAssay),
                          is.string(outputColumn))
  
  if (outputColumn %in% colnames(colData(object))) {
    warning("object already has a column named ", outputColumn, 
            " in colData: It will be overwritten!")
  }
  
  ## Calculate number of expressed loci
  res <- quantifyStrandwiseDivergentLoci(loci, object, inputAssay=inputAssay)
  m <- assay(res$'-', inputAssay)
  p <- assay(res$'+', inputAssay)
  expressed <- sapply(colnames(m), function(i) m[,i]+p[,i]>unexpressed)
  if (requirebidirectional)
    expressed <- expressed & sapply(colnames(m), 
                                    function(i) m[,i]>0 & p[,i]>0)
  colData(object)[, outputColumn] <- colSums(expressed)
  
  ## Return
  object
}