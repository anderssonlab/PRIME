require("CAGEfightR")
require("assertthat")

source("CAGEfightR_extensions/utils.R")

## object: SummarizedExperiment
## inputAssay: assay to normalize
## outputAssay: where to store normalized results
## conditionalColumn: numeric row-wise vector to normalize against
## offsetAssay: where to store offset

## GC normalization based on approach described in Pickrell et al 2010, Nature

conditionalNormalize <- function(object, inputAssay="counts", outputAssay="normalized", conditionalColumn="GC", offsetAssay=NULL, bins=200,
                                 sizeFactors=NULL, minCount=1, minSamples=1,aggregate.fn=sum) {
    
    assert_that(methods::is(object, "SummarizedExperiment"),
                inputAssay %in% assayNames(object),
                conditionalColumn %in% colnames(mcols(object)))

    message("binning data according to conditional column...")

    y <- assay(object,inputAssay)
    keep <- which(Matrix::rowSums(y>minCount) >= minSamples)

    x <- as.numeric(mcols(object)[,conditionalColumn])
    b <- as.character(Hmisc::cut2(x, unique(quantile(x, seq(1/bins, 1, 1/bins)))))
    b.m <- aggregate(x, by=list(b), FUN=mean)
    b.m.keep <- aggregate(x[keep], by=list(b[keep]), FUN=mean)
    n <- b.m[,1]
    b.m <- b.m[,2]
    names(b.m) <- n
    n <- b.m.keep[,1]
    b.m.keep <- b.m.keep[,2]
    names(b.m.keep) <- n

    y.b <- aggregate(as.matrix(y[keep,]),by=list(b[keep]),FUN=aggregate.fn)
    b.g <- y.b[,1]
    rownames(y.b) <- b.g
    y.b <- y.b[names(b.m.keep),]
    y.b <- data.matrix(y.b[,-1])

    message("calculating relative enrichments...")

    if (is.null(sizeFactors))
        sizeFactors <- rep(1,ncol(y))
    else
        sizeFactors <- sizeFactors / mean(sizeFactors)

    y.b <- sapply(1:ncol(y.b), function(i) y.b[,i] / sizeFactors[i])
    s <- colSums(y.b) / sum(y.b)
    
    f <- log2( t( t(y.b / rowSums(y.b)) / s ) )
    fit <- lapply(1:ncol(y), function(i) {
        b.m.i <- b.m.keep
        f.i <- f[,i]
        missing <- which(is.na(f[,i]) | is.infinite(f[,i]))
        if (length(missing)>0) {
            b.m.i <- b.m.i[-missing]
            f.i <- f.i[-missing]
        }
        smooth.spline(b.m.i, f.i, spar=1)
    })

    message("normalizing data according to enrichments...")

    offset <- -1*Matrix::Matrix(sapply(1:ncol(y), function(i) predict(fit[[i]],b.m[b])$y))
    y_normed <- Matrix::Matrix(round(y * 2^(offset)))
    y_normed <- Matrix::Matrix(sapply(1:ncol(y), function(i) {x <- y_normed[,i]; x[x==0 & y[,i]>0] <- 1; x}))

    dimnames(y_normed) <- dimnames(y)
    dimnames(offset) <- dimnames(y)
    
    assay(object, outputAssay) <- y_normed
    if (!is.null(offsetAssay))
        assay(object, offsetAssay) <- offset

    ## return
    object
}

normalizeBySizeFactors <- function(object, sizeFactors, inputAssay="counts", outputAssay="normalized") {

    assert_that(methods::is(object, "SummarizedExperiment"),
                inputAssay %in% assayNames(object),
                length(sizeFactors) == ncol(object))

    assay(object, outputAssay) <- Matrix::t(Matrix::t(assay(object, inputAssay)) / sizeFactors)

    object
}

TPMnormalizeBySizeFactors <- function(object, sizeFactors, inputAssay="counts", outputAssay="TPM") {

    assert_that(methods::is(object, "SummarizedExperiment"),
                inputAssay %in% assayNames(object),
                length(sizeFactors) == ncol(object))

    ## Centre size factors
    sf <- sizeFactors / mean(sizeFactors)

    ## Calculate scaling factors
    a <- assay(object, inputAssay)
    sf <- 1e6 / (sf * mean(Matrix::colSums(a)))
    
    assay(object, outputAssay) <- Matrix::t(Matrix::t(a) * sf)

    object
}
