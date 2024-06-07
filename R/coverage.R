require("CAGEfightR")
require("assertthat")

source("CAGEfightR_extensions/utils.R")

## ctss: RangedSummarizedExperiment
## categorical: factor

coverageAssociation <- function(ctss, categories, inputAssay="counts", outputColumn="CramersV", unexpressed=0) {
    
    ## dimensions
    assert_that(ncol(ctss)==length(categories),
                length(unique(categories))>1)

    if (outputColumn %in% colnames(rowData(ctss))) {
        warning("ctss already has a column named ", outputColumn, " in rowData: It will be overwritten!")
    }

    a <- as(assay(ctss,inputAssay)>unexpressed,"dgCMatrix")

    message("Calculating bias-corrected Cramer's V on ", nrow(a), " observations...")
    V <- apply(a,1,unbiasedCramersV,categories)
    
    rowData(ctss)[, outputColumn] <- V

    ctss
}

unbiasedCramersV <- function(x,y) {
    
    num.cat <- length(unique(y))
    num.val <- length(unique(x))
    
    if (num.val == 1 || num.cat == 1)
        return(0)
    
    gof <- suppressWarnings(chisq.test(x,y,correct=FALSE)$statistic)
    num.obs <- length(x)
    phi <- gof / num.obs
    
    phi.corr <- max(0, phi - (num.val - 1) * (num.cat - 1) / (num.obs - 1))
    if (phi.corr == 0)
        return(0)
    
    num.val.corr <- num.val - (num.val - 1)^2 / (num.obs - 1)
    num.cat.corr <- num.cat - (num.cat - 1)^2 / (num.obs - 1)
    
    sqrt(phi.corr / min(num.val.corr - 1, num.cat.corr - 1))
}
