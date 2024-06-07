require("CAGEfightR")
require("assertthat")

source("CAGEfightR_extensions/utils.R")

## object: RangedSummarizedExperiment
## categorical: factor

calcBatchSupport <- function(object, batch, inputAssay="counts", outputColumn="batchSupport", unexpressed=0, minSamples=1) {
    
    ## dimensions
    assert_that(ncol(object)==length(batch),
                length(unique(batch))>1,
                minSamples>=1)

    if (outputColumn %in% colnames(rowData(object))) {
        warning("object already has a column named ", outputColumn, " in rowData: It will be overwritten!")
    }

    message("Calculating support per batch...")
    a <- do.call("cbind", lapply(unique(batch), function(b) Matrix::rowSums(assay(object,inputAssay)[,batch==b,drop=FALSE]>unexpressed)))

    message("Calculating overall batch support...")
    rowData(object)[, outputColumn] <- as.integer(Matrix::rowSums(a>=minSamples))

    object
}
