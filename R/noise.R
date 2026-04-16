#' Estimate tag cluster noise in a given CTSS dataset by sampling and
#' quantifying regions of a fixed size across the genome.
#' 
#' @param object A \code{RangedSummarizedExperiment} object.
#' @param mask A \code{GRanges} object with the masked regions.
#' @param mappable A \code{GRanges} object with the mappable regions.
#' @param map_frac The fraction of bases that need to be mappable.
#' @param win_size The size of the windows to use.
#' @param num_win The number of windows to sample.
#' @param strand The strand to consider.
#' @param inputAssay The assay to use.
#' @param quantiles The quantiles to calculate.
#' 
#' @return A matrix with the quantiles for each sample.
#' 
#' @importFrom GenomicRanges tileGenome findOverlaps `strand<-` seqinfo
#' @importFrom CAGEfightR quantifyClusters
#' @importFrom GenomeInfoDb seqinfo seqlengths
#' @importFrom IRanges pintersect width
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom methods is
#' @importFrom assertthat assert_that is.string
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom stats quantile
#' 
#' @export
estimateNoise <- function(object, mask, mappable, map_frac=0.5, win_size=200, 
                          num_win=1e6, strand="+", inputAssay="counts", 
                          quantiles=c(0.9,0.95,0.99,0.999,0.9999,0.99999)) {
  
  assertthat::assert_that(
  methods::is(object, "RangedSummarizedExperiment"),
  methods::is(mask, "GRanges"),
  methods::is(mappable, "GRanges"),
  is.numeric(map_frac), length(map_frac) == 1, map_frac > 0, map_frac <= 1,
  is.numeric(win_size), length(win_size) == 1, win_size > 0,
  is.numeric(num_win), length(num_win) == 1, num_win > 0,
  assertthat::is.string(strand),
  strand %in% c("+", "-"),
  assertthat::is.string(inputAssay),
  inputAssay %in% SummarizedExperiment::assayNames(object),
  is.numeric(quantiles),
  all(quantiles > 0 & quantiles < 1)
)

  message("Creating unmasked windows...")
  ## Create tiling windows
  genome_gr <- GenomicRanges::tileGenome(GenomeInfoDb::seqinfo(object),tilewidth=win_size,
                          cut.last.tile.in.chrom=TRUE)
  ## Remove masked regions
  genome_gr <- IRanges::subsetByOverlaps(genome_gr, mask, maxgap=-1, type="any", 
                                invert=TRUE)
  
  message("Filtering windows based on mappability...")
  hits <- GenomicRanges::findOverlaps(mappable, genome_gr, maxgap=-1, type="any")
  olaps <- IRanges::pintersect(mappable[S4Vectors::queryHits(hits)], genome_gr[S4Vectors::subjectHits(hits)])
  df <- data.frame(subjectHits=S4Vectors::subjectHits(hits),coverage=IRanges::width(olaps) / 
                     IRanges::width(genome_gr[S4Vectors::subjectHits(hits)]))
  df <- aggregate(df, by=list(df$subjectHits), FUN=sum)
  
  ## Keep windows with at least map_frac fraction of bases uniquely mappable
  df <- subset(df, coverage >= map_frac)
  genome_gr <- genome_gr[df$Group.1]
  
  ## Sample windows randomly
  genome_gr <- genome_gr[sort(sample.int(length(genome_gr),
                                         min(num_win,length(genome_gr)),
                                         replace=FALSE))]
  
  message("Quantifying expression across samples...")
  GenomicRanges::strand(genome_gr) <- strand
  expr <- suppressMessages(SummarizedExperiment::assay(CAGEfightR::quantifyClusters(object, 
                                                  genome_gr, 
                                                  inputAssay = inputAssay),
                                 inputAssay))
  
  message("Calculating statistics...")
  res <- apply(expr,2,quantile,prob=quantiles)
  
  res
}

#' Estimate divergent loci noise in a given CTSS dataset by sampling and
#' quantifying regions of a fixed size across the genome.
#' 
#' @param object A \code{RangedSummarizedExperiment} object.
#' @param mask A \code{GRanges} object with the masked regions.
#' @param mappable_minus A \code{GRanges} object with the mappable regions for 
#' the minus strand.
#' @param mappable_plus A \code{GRanges} object with the mappable regions for 
#' the plus strand.
#' @param map_frac The fraction of bases that need to be mappable.
#' @param win_size The size of the windows to use.
#' @param num_win The number of windows to sample.
#' @param inputAssay The assay to use.
#' @param quantiles The quantiles to calculate.
#' 
#' @return A matrix with the quantiles for each sample.
#' 
#' @export
#' 
#' @importFrom GenomicRanges tileGenome findOverlaps `strand<-` seqinfo start end `start<-` `end<-`
#' @importFrom CAGEfightR quantifyClusters
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom IRanges pintersect width
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom methods is
#' @importFrom assertthat assert_that is.string
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom stats quantile
#' 
estimateDivergentNoise <- function(object, mask, mappable_minus, mappable_plus, 
                                   map_frac=0.5, win_size=200, num_win=1e6, 
                                   inputAssay="counts", 
                                   quantiles=c(0.9,0.95,0.99,0.999,0.9999,
                                               0.99999)) {
  
assertthat::assert_that(
  methods::is(object, "RangedSummarizedExperiment"),
  methods::is(mask, "GRanges"),
  methods::is(mappable_minus, "GRanges"),
  methods::is(mappable_plus, "GRanges"),
  is.numeric(map_frac), length(map_frac) == 1, map_frac > 0, map_frac <= 1,
  is.numeric(win_size), length(win_size) == 1, win_size > 0,
  is.numeric(num_win), length(num_win) == 1, num_win > 0,
  assertthat::is.string(inputAssay),
  inputAssay %in% SummarizedExperiment::assayNames(object),
  is.numeric(quantiles),
  all(quantiles > 0 & quantiles < 1)
)

  message("Creating unmasked mappable windows...")
  ## Create tiling windows
  genome_gr <- GenomicRanges::tileGenome(GenomeInfoDb::seqinfo(object),tilewidth=win_size*2+1,
                          cut.last.tile.in.chrom=TRUE)
  ## Remove masked regions
  genome_gr <- IRanges::subsetByOverlaps(genome_gr, mask, maxgap=-1, type="any", 
                                invert=TRUE)
  
  ## Create strand specific windows
  win_minus <- genome_gr
  win_plus <- genome_gr
  GenomicRanges::end(win_minus) <- GenomicRanges::start(win_minus) + win_size - 1
  GenomicRanges::start(win_plus) <- GenomicRanges::end(win_plus) - win_size + 1
  GenomicRanges::strand(win_minus) <- "-"
  GenomicRanges::strand(win_plus) <- "+"
  
  ## Require some mappability
  keep <- intersect(S4Vectors::subjectHits(GenomicRanges::findOverlaps(mappable_minus, win_minus, 
                                             maxgap=-1, type="any")),
                    S4Vectors::subjectHits(GenomicRanges::findOverlaps(mappable_plus, win_plus, 
                                             maxgap=-1, type="any")))
  win_minus <- win_minus[keep]
  win_plus <- win_plus[keep]
  
  message("Filtering windows based on mappability strand-wise...")
  hits_minus <- GenomicRanges::findOverlaps(mappable_minus, win_minus, maxgap=-1, type="any")
  hits_plus <- GenomicRanges::findOverlaps(mappable_plus, win_plus, maxgap=-1, type="any")
  
  olaps_minus <- IRanges::pintersect(mappable_minus[S4Vectors::queryHits(hits_minus)], 
                            win_minus[S4Vectors::subjectHits(hits_minus)])
  olaps_plus <- IRanges::pintersect(mappable_plus[S4Vectors::queryHits(hits_plus)], 
                           win_plus[S4Vectors::subjectHits(hits_plus)])
  
  df_minus <- data.frame(subjectHits=S4Vectors::subjectHits(hits_minus),
                         coverage=IRanges::width(olaps_minus) / 
                           IRanges::width(win_minus[S4Vectors::subjectHits(hits_minus)]))
  df_minus <- aggregate(df_minus, by=list(df_minus$subjectHits), FUN=sum)
  
  df_plus <- data.frame(subjectHits=S4Vectors::subjectHits(hits_plus),
                        coverage=IRanges::width(olaps_plus) / 
                          IRanges::width(win_plus[S4Vectors::subjectHits(hits_plus)]))
  df_plus <- aggregate(df_plus, by=list(df_plus$subjectHits), FUN=sum)
  
  ## Keep windows with at least map_frac fraction of bases uniquely mappable
  df_minus <- subset(df_minus, coverage >= map_frac)
  df_plus <- subset(df_plus, coverage >= map_frac)
  
  keep <- intersect(df_minus$Group.1,df_plus$Group.1)
  win_minus <- win_minus[keep]
  win_plus <- win_plus[keep]
  
  ## Sample windows randomly
  samp <- sample.int(length(keep),min(num_win,length(keep)),replace=FALSE)
  win_minus <- win_minus[samp]
  win_plus <- win_plus[samp]
  
  message("Quantifying expression strand_wise across samples...")
  expr_minus <- suppressMessages(SummarizedExperiment::assay(
    CAGEfightR::quantifyClusters(object, win_minus, inputAssay = inputAssay),inputAssay))
  expr_plus <- suppressMessages(SummarizedExperiment::assay(
    CAGEfightR::quantifyClusters(object, win_plus, inputAssay = inputAssay),inputAssay))
  
  expr <- expr_minus + expr_plus
  
  message("Calculating statistics...")
  res <- apply(expr,2,quantile,prob=quantiles)
  
  res
}
