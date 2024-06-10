#' Calculate the auto correlation of CTSSs over tag clusters according to 
#' pooled values and the cross correlation with opposite strand CTSSs
#'
#' @param loci A \code{GRanges} object of loci to focus on.
#' @param ctss A \code{RangedSummarizedExperiment} object containing the CTSSs.
#' @param max_dist The number of base pairs to consider around each locus.
#' 
#' @return A data frame containing the cumulative fraction of CTSSs around
#' the loci, for each sample.
#'
#' @export
cumulativeFractionAroundLoci <- function(loci, ctss, max_dist=1000) {
  
  p <- GRanges(seqnames=seqnames(loci),loci$thick,strand="+")
  m <- GRanges(seqnames=seqnames(loci),loci$thick,strand="-")
  t <- as(rowRanges(ctss),"GRanges")
  
  message("findining nearest loci strand-wise...")
  p_nearest <- follow(t, p, select="last")
  m_nearest <- follow(t, m, select="last")
  
  message("calculating distances...")
  d <- rep(Inf,length(t))
  non.na <- which(!is.na(p_nearest))
  d[non.na] <- distance(t[non.na],p[p_nearest[non.na]])
  non.na <- which(!is.na(m_nearest))
  d[non.na] <- distance(t[non.na],m[m_nearest[non.na]])
  
  if (!"totalTags" %in% colnames(colData(ctss))) {
    message("calculating number of tags per sample...")
    ctss <- CAGEfightR::calcTotalTags(ctss)
  }
  
  message("calculating cumulative fractions...")
  focus <- which(d<=max_dist)
  d <- d[focus]
  ctss <- ctss[focus]
  
  o <- order(d)
  d <- d[o]
  fracs <- sapply(colnames(ctss), function(n) {
    cat(".")
    a <- as.numeric(assay(ctss,"counts")[o,n])
    nz <- which(a>0)
    da <- d[nz]
    a <- a[nz]
    sapply(seq(0,max_dist,by=1), function(i) sum(a[da==i]))
  })
  cat("\n")
  fracs <- sapply(colnames(fracs),function(n) {
    cumsum(fracs[,n]/ctss$totalTags[n]))
  })

cumfrac <- cbind(data.frame(distance=seq(0,max_dist,by=1),fracs))

cumfrac
}
