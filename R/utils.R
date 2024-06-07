
calcNumberCTSSs <- function(object, inputAssay = "counts", outputColumn = "numberCTSSs", unexpressed = 0) {
    ## Prechecks
    assert_that(methods::is(object, "SummarizedExperiment"),
                not_empty(object),
                inputAssay %in% assayNames(object),
                is.string(inputAssay),
                is.string(outputColumn))

    if (outputColumn %in% colnames(colData(object))) {
        warning("object already has a column named ", outputColumn, " in colData: It will be overwritten!")
    }

    ## Calculate number of unique CTSS positions
    a <- assay(object, inputAssay)
    colData(object)[, outputColumn] <- sapply(colnames(a), function(i) sum(a[,i]>unexpressed))

    ## Return
    object
}

calcNumberGenes <- function(object, txModels, inputAssay = "counts", outputColumn = "numberGenes", unexpressed = 0) {
    ## Prechecks
    assert_that(methods::is(object, "SummarizedExperiment"),
                not_empty(object),
                inputAssay %in% assayNames(object),
                is.string(inputAssay),
                is.string(outputColumn))

    if (outputColumn %in% colnames(colData(object))) {
        warning("object already has a column named ", outputColumn, " in colData: It will be overwritten!")
    }

    if (!"geneID" %in% colnames(mcols(object)))
        object <- assignGeneID(object, geneModels = txModels, outputColumn = "geneID")
    genelevel <- quantifyGenes(object, genes="geneID", inputAssay=inputAssay)

    ## Calculate number of expressed genes
    a <- assay(genelevel, inputAssay)
    colData(object)[, outputColumn] <- sapply(colnames(a), function(i) sum(a[,i]>unexpressed))

    ## Return
    object
}

calcNumberDivergentLoci <- function(object, loci, inputAssay="counts", outputColumn = "numberLoci", unexpressed = 0, requirebidirectional=FALSE) {
    ## Prechecks
    assert_that(methods::is(object, "SummarizedExperiment"),
                not_empty(object),
                inputAssay %in% assayNames(object),
                is.string(inputAssay),
                is.string(outputColumn))

    if (outputColumn %in% colnames(colData(object))) {
        warning("object already has a column named ", outputColumn, " in colData: It will be overwritten!")
    }

    ## Calculate number of expressed loci
    res <- quantifyStrandwiseDivergentLoci(loci, object, inputAssay=inputAssay)
    m <- assay(res$'-', inputAssay)
    p <- assay(res$'+', inputAssay)
    expressed <- sapply(colnames(m), function(i) m[,i]+p[,i]>unexpressed)
    if (requirebidirectional)
        expressed <- expressed & sapply(colnames(m), function(i) m[,i]>0 & p[,i]>0)
    colData(object)[, outputColumn] <- colSums(expressed)

    ## Return
    object
}

quantifyClustersOlap <- function(object, clusters, inputAssay) {

    revmap <- GenomicRanges::reduce(clusters, with.revmap=TRUE)$revmap
    max.nr <- max(sapply(revmap, length))
    ids <- lapply(1:max.nr, function(level) setdiff(sapply(revmap, `[`, level),NA))

    res <- BiocParallel::bplapply(1:max.nr, function(i) {
        clu <- clusters[ids[[i]]]
        obj <- object
        seqinfo(obj) <- seqinfo(clu)

        ##special case, CAGEfightR::quantifyClusters does not deal with granges objects of length 1
        spec <- FALSE
        if (length(clu) == 1) {
            spec <- TRUE
            clu <- c(clu, shift(clu,width(clu)))
        }

        suppressMessages(m_ <- CAGEfightR::quantifyClusters(
                                                object = obj,
                                                clusters = clu,
                                                inputAssay = inputAssay
                                            ))
        if (spec)
            m_ <- m_[1]
        SummarizedExperiment::assays(m_)[[inputAssay]]
    })

    mat <- do.call("rbind", res)[order(unlist(ids)),]

    rownames(mat) <- names(clusters)

    o <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(mat),
                                                    rowRanges = clusters,
                                                    colData = SummarizedExperiment::colData(object))

    SummarizedExperiment::assayNames(o) <- inputAssay

    o
}

## loci: GRanges
## ctss: RangedSummarisedExperiment
quantifyStrandwiseDivergentLoci <- function(loci, ctss, inputAssay = "counts", requireDisjoint=TRUE) {
  win_1 <- loci
  BiocGenerics::end(win_1) <- start(loci$thick) - 1
  BiocGenerics::strand(win_1) <- "-"

  win_2 <- loci
  BiocGenerics::start(win_2) <- start(loci$thick) + 1
  BiocGenerics::strand(win_2) <- "+"

  m1 <- matrix()
  if (!requireDisjoint && !isDisjoint(win_1)) {
      m1_ <- quantifyClustersOlap(
          object = ctss,
          clusters = win_1,
          inputAssay = inputAssay
      )
      m1 <- SummarizedExperiment::assays(m1_)[[inputAssay]]
  }
  else {
      m1_ <- CAGEfightR::quantifyClusters(
                             object = ctss,
                             clusters = win_1,
                             inputAssay = inputAssay
                         )
      m1 <- SummarizedExperiment::assays(m1_)[[inputAssay]]
  }

  m2 <- matrix()
  if (!requireDisjoint && !isDisjoint(win_2)) {
      m2_ <- quantifyClustersOlap(
          object = ctss,
          clusters = win_2,
          inputAssay = inputAssay
      )
      m2 <- SummarizedExperiment::assays(m1_)[[inputAssay]]
  }
  else {
      m2_ <- CAGEfightR::quantifyClusters(
                             object = ctss,
                             clusters = win_2,
                             inputAssay = inputAssay
                         )
      m2 <- SummarizedExperiment::assays(m2_)[[inputAssay]]
  }

  m <- SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(m1),
    rowRanges = loci,
    colData = SummarizedExperiment::colData(ctss)
  )
  SummarizedExperiment::assayNames(m) <- inputAssay
  p <- SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(m2),
    rowRanges = loci,
    colData = SummarizedExperiment::colData(ctss)
  )
  SummarizedExperiment::assayNames(p) <- inputAssay

  res <- base::list(m,p)
  base::names(res) <- c("-","+")

  res
}

## loci: GRanges
## ctss: RangedSummarisedExperiment
quantifyDivergentLoci <- function(loci, ctss, inputAssay="counts", requireDisjoint=TRUE) {

  res <- quantifyStrandwiseDivergentLoci(loci, ctss, inputAssay, requireDisjoint)

  o <- SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(
      SummarizedExperiment::assays(res$'-')[[inputAssay]] +
        SummarizedExperiment::assays(res$'+')[[inputAssay]]
    ),
    rowRanges = loci,
    colData = SummarizedExperiment::colData(ctss)
  )
  SummarizedExperiment::assayNames(o) <- inputAssay

  o
}

## nonzero function from the DAPAR package (https://rdrr.io/bioc/DAPAR/src/R/utils.R)
nonzero <- function(x){
    stopifnot(inherits(x, "dgCMatrix"))
    if (all(x@p == 0))
        return(matrix(0, nrow=0, ncol=2,
                      dimnames=list(character(0), c("row","col"))))
    res <- cbind(x@i+1, rep(seq(dim(x)[2]), diff(x@p)))
    colnames(res) <- c("row", "col")
    res <- res[x@x != 0, , drop = FALSE]
    return(res)
}

### Helper functions not exported by CAGEfightR

TCstats <- function(coverage_plus, coverage_minus, tcs_plus, tcs_minus) {
                                        # Check classes
    stopifnot(methods::is(coverage_plus, "SimpleRleList"),
              methods::is(coverage_minus, "SimpleRleList"),
              methods::is(tcs_plus, "CompressedIRangesList"),
              methods::is(tcs_minus, "CompressedIRangesList"))

                                        # Check seqlevels
    stopifnot(length(unique(list(names(coverage_plus),
                                 names(tcs_plus),
                                 names(coverage_minus),
                                 names(tcs_minus)))) == 1)

                                        # Obtain views
    views_plus <- Views(coverage_plus, tcs_plus)
    views_minus <- Views(coverage_minus, tcs_minus)

                                        # Calculate Sums
    sum_plus <- unlist(viewSums(views_plus))
    sum_minus <- unlist(viewSums(views_minus))

                                        # Find peaks
    ranges_plus <- viewRangeMaxs(views_plus)
    ranges_minus <- viewRangeMaxs(views_minus)
    ranges_plus <- resize(unlist(ranges_plus), width = 1, fix = "center")
    ranges_minus <- resize(unlist(ranges_minus), width = 1, fix = "center")

                                        # Merge into GRanges
    TCs <- c(GRanges(tcs_plus,
                     strand = "+",
                     score = sum_plus,
                     thick = ranges_plus),
             GRanges(tcs_minus,
                     strand = "-",
                     score = sum_minus,
                     thick = ranges_minus))

                                        # Names as IDs for both ranges and peaks
    TC_ids <- paste0(seqnames(TCs), ":", start(TCs), "-", end(TCs), ";", strand(TCs))
    names(TCs) <- TC_ids
    names(TCs$thick) <- TC_ids

                                        # Return
    TCs
}

summarizeWidths <- function(gr) {
                                        # Checks
    stopifnot(methods::is(gr, "GRanges"))

                                        # Cut up widths
    x <- cut(width(gr), breaks = c(1, 10, 100, 1000, Inf), labels = c(">=1", ">=10",
                                                                      ">=100", ">=1000"), include.lowest = TRUE)

                                        # Get freqs and props
    y <- table(Width = x)
    z <- prop.table(y)

                                        # Format to data.frame
    w <- merge(as.data.frame(y, responseName = "Count"),
               as.data.frame(z, responseName = "Percent"))

                                        # Add Total row
    w <- rbind(data.frame(Width = "Total",
                          Count = sum(w$Count),
                          Percent = sum(w$Percent)), w)

                                        # Reformat to percent
    w$Percent <- paste0(format(w$Percent * 100, digits = 1), " %")

                                        # To string and message
    s <- paste(utils::capture.output(print(w, row.names = FALSE)), collapse = "\n")
    message(s)
    }

splitByStrand <- function(object) {
    split(object, strand(object))
}

splitPooled <- function(object, weight="score"){

    ## Split by strand
    o <- splitByStrand(object)

    ## Calculate coverage
    o <- lapply(o, coverage, weight=weight)

    ## Round to handle floating point errors
    o <- lapply(o, round, digits=9)

    ## Return
    o
}
