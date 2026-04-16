#' Calculate the cumulative fraction of CTSSs around loci.
#'
#' @param loci A \code{GRanges} object of loci to focus on.
#' @param ctss A \code{RangedSummarizedExperiment} object containing the CTSSs.
#' @param max_dist The number of base pairs to consider around each locus.
#'
#' @return A data frame containing the cumulative fraction of CTSSs around
#' the loci, for each sample.
#'
#' @export
#'
#' @import GenomicRanges
#' @import CAGEfightR
#' @importFrom IRanges distance follow
#' @importFrom SummarizedExperiment assay colData rowRanges
#' @importFrom assertthat assert_that is.numeric is.flag
#' @importFrom S4Vectors mcols
#'
#'
cumulativeFractionAroundLoci <- function(loci, ctss, max_dist = 1000) {
  
  assertthat::assert_that(
  methods::is(loci, "GRanges"),
  "thick" %in% colnames(S4Vectors::mcols(loci)),
  methods::is(ctss, "RangedSummarizedExperiment"),
  is.numeric(max_dist), length(max_dist) == 1, max_dist >= 0,
  "counts" %in% SummarizedExperiment::assayNames(ctss) 
)

  p <- GRanges(seqnames = seqnames(loci), loci$thick, strand = "+")
  m <- GRanges(seqnames = seqnames(loci), loci$thick, strand = "-")
  t <- as(rowRanges(ctss), "GRanges")

  message("findining nearest loci strand-wise...")
  p_nearest <- follow(t, p, select = "last")
  m_nearest <- follow(t, m, select = "last")

  message("calculating distances...")
  d <- rep(Inf, length(t))
  non.na <- which(!is.na(p_nearest))
  d[non.na] <- distance(t[non.na], p[p_nearest[non.na]])
  non.na <- which(!is.na(m_nearest))
  d[non.na] <- distance(t[non.na], m[m_nearest[non.na]])

  if (!"totalTags" %in% colnames(colData(ctss))) {
    message("calculating number of tags per sample...")
    ctss <- calcTotalTags(ctss)
  }

  message("calculating cumulative fractions...")
  focus <- which(d <= max_dist)
  d <- d[focus]
  ctss <- ctss[focus]

  o <- order(d)
  d <- d[o]
  fracs <- sapply(colnames(ctss), function(n) {
    cat(".")
    a <- as.numeric(assay(ctss, "counts")[o, n])
    nz <- which(a > 0)
    da <- d[nz]
    a <- a[nz]
    sapply(seq(0, max_dist, by = 1), function(i) sum(a[da == i]))
  })
  cat("\n")
  fracs <- sapply(colnames(fracs), function(n) {
    cumsum(fracs[, n] / ctss$totalTags[n])
  })

  cumfrac <- cbind(data.frame(distance = seq(0, max_dist, by = 1), fracs))

  cumfrac
}

#' Calculate the cumulative fraction of CTSSs around loci, considering unique
#' tags and the fraction of loci with CTSSs.
#'
#' @param loci A \code{GRanges} object of loci to focus on.
#' @param ctss A \code{RangedSummarizedExperiment} object containing the CTSSs.
#' @param max_dist The number of base pairs to consider around each locus.
#' @param inputAssay The assay in the \code{ctss} object to use for counting tags (default is "counts").
#'
#' @return A list of three data frames: cumulative fraction of tags, cumulative fraction of unique tags, and cumulative fraction of loci with CTSSs, for each sample.
#' @export
#' @import GenomicRanges
#' @import CAGEfightR
#' @importFrom IRanges distance follow
#' @importFrom SummarizedExperiment assay colData rowRanges
#' @importFrom assertthat assert_that is.numeric is.flag is.string
#' @importFrom S4Vectors mcols
#' 
#' 
cumulativeFractionsAroundLoci <- function(loci, ctss, max_dist = 1000, inputAssay = "counts") {
   
  assertthat::assert_that(
  methods::is(loci, "GRanges"),
  "thick" %in% colnames(S4Vectors::mcols(loci)),
  methods::is(ctss, "RangedSummarizedExperiment"),
  is.numeric(max_dist), length(max_dist) == 1, max_dist >= 0,
  assertthat::is.string(inputAssay),
  inputAssay %in% RangedSummarizedExperiment::assayNames(ctss)
)

  g <- GenomicRanges::GRanges(
    seqnames = GenomeInfoDb::seqnames(loci),
    ranges = loci$thick,
    strand = "*"
  )
  t <- methods::as(SummarizedExperiment::rowRanges(ctss), "GRanges")

  base::message("finding nearest loci...")
  g_nearest <- IRanges::nearest(t, g)

  base::message("calculating distances...")
  d <- base::rep(Inf, base::length(t))
  non.na <- base::which(!base::is.na(g_nearest))
  d[non.na] <- IRanges::distance(t[non.na], g[g_nearest[non.na]])

  base::message("calculating totalTags and unique Tags...")
  ctss <- base::suppressWarnings(
    CAGEfightR::calcTotalTags(
      object = ctss,
      inputAssay = "counts"
    )
  )
  SummarizedExperiment::colData(ctss)$uniqTags <- Matrix::colSums(
    SummarizedExperiment::assay(ctss, inputAssay = "counts") > 0
  )

  base::message("calculating cumulative fractions...")
  focus <- base::which(d <= max_dist)
  d <- d[focus]
  ctss <- ctss[focus]

  o <- base::order(d)
  d <- d[o]

  unique_loc_frac <- base::sapply(
    X = base::colnames(ctss),
    FUN = function(n) {
      base::cat(".")
      a <- base::as.numeric(
        SummarizedExperiment::assay(ctss, inputAssay)[o, n]
      )
      a[a > 0] <- 1
      nz <- base::which(a > 0)
      da <- d[nz]
      a_sub <- a[nz]
      return(
        base::sapply(
          X = base::seq(0, max_dist, by = 1),
          FUN = function(i) {
            base::sum(a_sub[da == i])
          }
        )
      )
    }
  )
  base::cat("\n")

  fracs <- base::sapply(
    X = base::colnames(ctss),
    FUN = function(n) {
      base::cat(".")
      a <- base::as.numeric(
        SummarizedExperiment::assay(ctss, inputAssay)[o, n]
      )
      nz <- base::which(a > 0)
      da <- d[nz]
      a_sub <- a[nz]
      return(
        base::sapply(
          X = base::seq(0, max_dist, by = 1),
          FUN = function(i) {
            base::sum(a_sub[da == i])
          }
        )
      )
    }
  )
  base::cat("\n")

  t <- methods::as(SummarizedExperiment::rowRanges(ctss), "GRanges")

  loc_fracs <- base::sapply(
    X = base::colnames(ctss),
    FUN = function(n) {
      base::cat(".")
      a <- base::as.numeric(
        SummarizedExperiment::assay(ctss, inputAssay)[o, n]
      )
      nz <- t[base::which(a > 0)]
      dist <- IRanges::distanceToNearest(g, nz)
      dist <- base::as.data.frame(dist)[, 3]
      return(
        base::sapply(
          X = base::seq(0, max_dist, by = 1),
          FUN = function(i) {
            base::sum(dist == i)
          }
        )
      )
    }
  )

  unique_loc_frac <- base::sapply(
    X = base::colnames(unique_loc_frac),
    FUN = function(n) {
      base::cumsum(unique_loc_frac[, n] / ctss$uniqTags[n])
    }
  )
  fracs <- base::sapply(
    X = base::colnames(fracs),
    FUN = function(n) {
      base::cumsum(fracs[, n] / ctss$totalTags[n])
    }
  )
  loc_fracs <- base::sapply(
    X = base::colnames(loc_fracs),
    FUN = function(n) {
      base::cumsum(loc_fracs[, n] / base::length(g))
    }
  )

  ucumfrac <- base::cbind(
    base::data.frame(
      distance = base::seq(0, max_dist, by = 1),
      unique_loc_frac
    )
  )
  cumfrac <- base::cbind(
    base::data.frame(
      distance = base::seq(0, max_dist, by = 1),
      fracs
    )
  )
  locifrac <- base::cbind(
    base::data.frame(
      distance = base::seq(0, max_dist, by = 1),
      loc_fracs
    )
  )

  fraclist <- base::list(cumfrac, ucumfrac, locifrac)
  return(fraclist)
}
