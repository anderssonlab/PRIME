require("CAGEfightR")

source("CAGEfightR_extensions/utils.R")

calcComplexity <- function(object, txModels, step=1e6, CTSSunexpressed=1, geneunexpressed=9, minCTSSsupport=2) {

    object <- suppressWarnings(calcTotalTags(object, inputAssay="counts"))
    object <- suppressMessages(suppressWarnings(assignGeneID(object, geneModels = txModels, outputColumn = "geneID")))

    targets <- seq(step,max(object$totalTags),by=step)

    message("Subsampling counts")
    res <- c(list(data.frame(target=0,sample=colnames(object),totalTags=0,numberCTSSs=0,numberGenes=0)),
             lapply(targets, function(t) {
                 message(t)

                 x <- subsampleTarget(object, "counts", t)
                 x <- suppressWarnings(calcTotalTags(x, inputAssay="counts"))
                 if (minCTSSsupport > 1)
                     x <- suppressMessages(subsetBySupport(x, unexpressed=CTSSunexpressed, minSamples=minCTSSsupport))

                 x <- suppressWarnings(calcNumberCTSSs(x, inputAssay="counts", unexpressed = CTSSunexpressed))
                 x <- suppressWarnings(calcNumberGenes(x, txModels, inputAssay="counts", unexpressed = geneunexpressed))

                 df <- data.frame(target=t,sample=colnames(x),totalTags=x$totalTags,numberCTSSs=x$numberCTSSs,numberGenes=x$numberGenes)
                 df[object$totalTags<t,c("totalTags","numberCTSSs","numberGenes")] <- NA

                 df
             }))

    do.call("rbind",res)
}

calcCTSSComplexity <- function(object, step=1e6, CTSSunexpressed=1, minCTSSsupport=2) {

    object <- suppressWarnings(calcTotalTags(object, inputAssay="counts"))

    targets <- seq(step,max(object$totalTags),by=step)

    message("Subsampling counts")
    res <- c(list(data.frame(target=0,sample=colnames(object),totalTags=0,numberCTSSs=0)),
             lapply(targets, function(t) {
                 message(t)

                 x <- subsampleTarget(object, "counts", t)
                 x <- suppressWarnings(calcTotalTags(x, inputAssay="counts"))
                 if (minCTSSsupport > 1)
                     x <- suppressMessages(subsetBySupport(x, unexpressed=CTSSunexpressed, minSamples=minCTSSsupport))

                 x <- suppressWarnings(calcNumberCTSSs(x, inputAssay="counts", unexpressed = CTSSunexpressed))

                 df <- data.frame(target=t,sample=colnames(x),totalTags=x$totalTags,numberCTSSs=x$numberCTSSs)
                 df[object$totalTags<t,c("totalTags","numberCTSSs")] <- NA

                 df
             }))

    do.call("rbind",res)
}

calcGeneComplexity <- function(object, txModels, step=1e6, CTSSunexpressed=1, geneunexpressed=9, minCTSSsupport=2) {

    object <- suppressWarnings(calcTotalTags(object, inputAssay="counts"))
    object <- suppressMessages(suppressWarnings(assignGeneID(object, geneModels = txModels, outputColumn = "geneID")))

    targets <- seq(step,max(object$totalTags),by=step)

    message("Subsampling counts")
    res <- c(list(data.frame(target=0,sample=colnames(object),totalTags=0,numberGenes=0)),
             lapply(targets, function(t) {
                 message(t)

                 x <- subsampleTarget(object, "counts", t)
                 x <- suppressWarnings(calcTotalTags(x, inputAssay="counts"))
                 if (minCTSSsupport > 1)
                     x <- suppressMessages(subsetBySupport(x, unexpressed=CTSSunexpressed, minSamples=minCTSSsupport))

                 x <- suppressWarnings(calcNumberGenes(x, txModels, inputAssay="counts", unexpressed = geneunexpressed))

                 df <- data.frame(target=t,sample=colnames(x),totalTags=x$totalTags,numberGenes=x$numberGenes)
                 df[object$totalTags<t,c("totalTags","numberGenes")] <- NA

                 df
             }))

    do.call("rbind",res)
}

calcDivergentLociComplexity <- function(object, loci, step=1e6, CTSSunexpressed=1, lociunexpressed=2, requirebidirectional=FALSE, minCTSSsupport=2) {

    object <- suppressWarnings(calcTotalTags(object, inputAssay="counts"))

    targets <- seq(step,max(object$totalTags),by=step)

    message("Subsampling counts")
    res <- c(list(data.frame(target=0,sample=colnames(object),totalTags=0,numberLoci=0)),
             lapply(targets, function(t) {
                 message(t)

                 x <- subsampleTarget(object, "counts", t)
                 x <- suppressWarnings(calcTotalTags(x, inputAssay="counts"))
                 if (minCTSSsupport > 1)
                     x <- suppressMessages(subsetBySupport(x, unexpressed=CTSSunexpressed, minSamples=minCTSSsupport))

                 x <- suppressWarnings(calcNumberDivergentLoci(x, loci=loci, inputAssay="counts", unexpressed=lociunexpressed, requirebidirectional=requirebidirectional))

                 df <- data.frame(target=t,sample=colnames(x),totalTags=x$totalTags,numberLoci=x$numberLoci)
                 df[object$totalTags<t,c("totalTags","numberLoci")] <- NA

                 df
             }))

    do.call("rbind",res)
}
