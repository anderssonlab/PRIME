require("CAGEfightR")
require("DescTools")

## Object: GRanges
## ctss: RangedSummarizedExperiment
CTSSdispersion <- function(object, ctss, inputAssay="TPM", outputColumn="TSS_GINI", fn, ...) {

    assertthat::assert_that(isSorted(ctss))
    assertthat::assert_that(isSorted(object))

    ## Find overlaps
    message("Finding overlaps...")

    hits <- findOverlaps(query = ctss,
                         subject = object,
                         select = "arbitrary")
    hits <- factor(hits, levels=seq_along(object))
    ids <- as.numeric(unlist(split(1:length(hits), hits)))
    ctss <- ctss[ids,]
    ids <- split(1:nrow(ctss), hits[-which(is.na(hits))])

    ## Calculate dispersion
    message("Calculating TSS usage dispersion...")

    data <- assay(ctss, inputAssay)
    value <- unlist(bplapply(ids, function(x) {
        fn(data[x,,drop=FALSE],...)
    }))

    mcols(object)[, outputColumn] <- value

    object
}

## object: RangedSummarizedExperiment
TCdispersion <- function(object, inputAssay="TPM", outputColumn="TC_GINI", fn=tc_gini, ...) {

    message("Calculating TC dispersion...")

    data <- assay(object, inputAssay)
    value <- unlist(bplapply(1:nrow(data), function(i) {
        fn(data[i,],...)
    }))

    mcols(object)[, outputColumn] <- value

    object
}

## object: RangedSummarizedExperiment
calcMedian <- function(object, inputAssay="counts", outputColumn="median") {

    data <- assay(object, inputAssay)
    value <- apply(data,1,median,na.rm=TRUE)

    mcols(object)[, outputColumn] <- value

    object
}

## object: RangedSummarizedExperiment
calcMax <- function(object, inputAssay="counts", outputColumn="median") {

    data <- assay(object, inputAssay)
    value <- apply(data,1,max,na.rm=TRUE)

    mcols(object)[, outputColumn] <- value

    object
}

## object: RangedSummarizedExperiment
calcMean <- function(object, inputAssay="counts", outputColumn="mean") {

    data <- assay(object, inputAssay)
    value <- Matrix::rowMeans(data,na.rm=TRUE)

    mcols(object)[, outputColumn] <- value

    object
}

TCdispersionData <- function(object, ctss, inputAssay="TPM", fn=tss_disp_data, ...) {

    assertthat::assert_that(isSorted(ctss))
    assertthat::assert_that(isSorted(object))

    ## Find overlaps
    message("Finding overlaps...")

    hits <- findOverlaps(query = ctss,
                         subject = object,
                         select = "arbitrary")
    hits <- factor(hits, levels=seq_along(object))
    ids <- as.numeric(unlist(split(1:length(hits), hits)))
    ctss <- ctss[ids,]
    ids <- split(1:nrow(ctss), hits[-which(is.na(hits))])

    ## Calculate dispersion
    message("Extracting data...")

    data <- assay(ctss, inputAssay)
    value <- bplapply(ids, function(x) {
        fn(data[x,,drop=FALSE],...)
    })

    value
}

tss_disp_data <- function(d, pseudo=rep(0,ncol(d)), max_zero_prop=1) {
    d <- as.matrix(d)
    if (nrow(d)==1)
        return(0)
    if (max_zero_prop < 1)
    {
        propz <- apply(d,1,function(x) sum(x==0))/ncol(d)
        idx <- which(propz <= max_zero_prop)
    }
    d.pseudo <- t(t(d)+pseudo)
    dnorm <- apply(d.pseudo,2,function(x) x/sum(x))
    list(expr=colSums(d)+pseudo,
         orig=d.pseudo[idx,,drop=FALSE],
         norm=dnorm[idx,,drop=FALSE])
}

rowVars <- function(x, ...) {
    sqr <- function(x)  x*x
    n <- Matrix::rowSums(!is.na(x))
    n[n<=1] <- NA
    return(rowSums(sqr(x-Matrix::rowMeans(x, ...)), ...)/(n-1))
}

rowSds <- function(x, ...)
      sqrt(rowVars(x, ...))

calcCV <- function(object, inputAssay="counts", outputColumn="CV") {

    data <- assay(object, inputAssay)
    value <- rowSds(data,na.rm=TRUE)/Matrix::rowMeans(data,na.rm=TRUE)

    mcols(object)[, outputColumn] <- value

    object
}

associateWithTC <- function(object, ctss) {

    assertthat::assert_that(isSorted(ctss))
    assertthat::assert_that(isSorted(object))

    hits <- findOverlaps(query = ctss,
                         subject = object,
                         select = "arbitrary")
    hits <- factor(hits, levels=seq_along(object))
    ids <- as.numeric(unlist(split(1:length(hits), hits)))
    ctss <-  ctss[ids,]
    mcols(ctss)[,"TC"] <- rownames(object)[hits[-which(is.na(hits))]]

    ctss
}

DMadjustedCV <- function(object, inputAssay="TPM", prefix="") {

    require(zoo)

    object <- calcCV(object, inputAssay, paste(prefix,"CV",sep=""))
    object <- calcMean(object, inputAssay, paste(prefix,"mean",sep=""))
    mcols(object)[,paste(prefix,"log10_CV2",sep="")] <- log10(rowData(object)[,paste(prefix,"CV",sep="")]^2)

    mean_order <- order(rowData(object)[,paste(prefix,"mean",sep="")])

    rolling_median <- rollapply(mcols(object)[,paste(prefix,"log10_CV2",sep="")][mean_order],width=50,by=25,FUN=median,fill=list("extend","extend","NA"))
    nas <- which(is.na(rolling_median))
    rolling_median[nas] <- median(mcols(object)[,paste(prefix,"log10_CV2",sep="")][mean_order][nas])
    names(rolling_median) <- (1:length(object))[mean_order]
    reorder <- match(1:length(object),names(rolling_median))
    rolling_median <- rolling_median[reorder]
    mcols(object)[, paste(prefix,"roll_median_log10_CV2",sep="")] <- rolling_median

    mcols(object)[, paste(prefix,"adjusted_log10_CV2",sep="")] <- mcols(object)[,paste(prefix,"log10_CV2",sep="")] - rolling_median

    object
}

quantifyTCfraction <- function(ctss, TCcolumn="TC", inputAssay="TPM", outputAssay="fraction", pseudo=1e6/ctss$totalTags) {

    assertthat::assert_that(isSorted(ctss))

    data <- assay(ctss,inputAssay)
    TCs <- mcols(ctss)[,TCcolumn]

    expr <- Matrix::fac2sparse(factor(TCs,levels=unique(TCs)))
    expr <- expr %*% data
    expr <- expr + pseudo
    expr <- expr[TCs,]

    assay(ctss, outputAssay) <- data / expr

    ctss
}

## Object: GRanges
## ctss: RangedSummarizedExperiment
numCTSSs <- function(object, ctss, inputAssay="TPM", outputColumn="num_CTSSs", max_zero_prop=1) {

    assertthat::assert_that(isSorted(ctss))
    assertthat::assert_that(isSorted(object))

    ## Find overlaps
    message("Finding overlaps...")

    hits <- findOverlaps(query = ctss,
                         subject = object,
                         select = "arbitrary")
    hits <- factor(hits, levels=seq_along(object))
    ids <- as.numeric(unlist(split(1:length(hits), hits)))
    ctss <- ctss[ids,]
    ids <- split(1:nrow(ctss), hits[-which(is.na(hits))])

    ## Calculate dispersion
    message("Calculating number of CTSSs with valid zero proportion...")

    data <- assay(ctss, inputAssay)
    value <- unlist(bplapply(ids, function(x) {
        x <- as.matrix(data[x,,drop=FALSE])
        sum(apply(x,1,function(y) sum(y==0))/ncol(x) <= max_zero_prop)
    }))

    mcols(object)[, outputColumn] <- value

    object
}

tc_gini <- function(d, pseudo=rep(0,length(d)), max_zero_prop=0.9) {
    d <- as.numeric(d)
    if (max_zero_prop < 1)
    {
        propz <- sum(d==0)/length(d)
        if (propz > max_zero_prop)
            return(NA)
    }
    d <- d+pseudo
    Gini(d)
}

madm <- function(x) {
    m <- median(x,na.rm=TRUE)
    if (m==0)
        return(NA)
    median(abs(x-m),na.rm=TRUE)/m
}

tc_madm <- function(d, pseudo=rep(0,length(d)), max_zero_prop=0.5) {
    d <- as.numeric(d)
    if (max_zero_prop < 1)
    {
        propz <- sum(d==0)/length(d)
        if (propz > max_zero_prop)
            return(NA)
    }
    d <- d+pseudo
    madm(d)
}
