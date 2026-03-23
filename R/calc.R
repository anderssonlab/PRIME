
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

calcAnnoCTSS <- function(data, txModels, uniqueCTSS = FALSE) {
    
    annos <- CAGEfightR::assignTxType(
        object = data,
        txModels = txModels
      )
    
    #Identify all txTypes
    types <- base::unique(SummarizedExperiment::rowData(annos)$txType)
    
    # Initialize list to store total tags for each txType
    totalTags_list <- base::vector("list", length = base::length(types))
    
    # Loop through each txType
    for (i in base::seq_along(types)) {
      
        # Subset data based on current txType
        current_subset <- S4Vectors::subset(annos, txType %in% types[i])
        if(uniqueCTSS){
          
            # Calculate totalTags for current subset
            current_subset <- Matrix::colSums(
                SummarizedExperiment::assay(current_subset) > 0
              ) 
            totalTags_list[[i]] <- current_subset
            
          } else {
            
            # Calculate totalTags for current subset
            current_subset <- SummarizedExperiment::colData(
                CAGEfightR::calcTotalTags(current_subset)
              )
            totalTags_list[[i]] <- current_subset$totalTags
            
          }
      }
    
    # Combine totalTags into a data frame
    ann_counts <- base::as.data.frame(
        base::do.call(base::cbind, totalTags_list)
      )
    base::colnames(ann_counts) <- types
    return(ann_counts)
  }