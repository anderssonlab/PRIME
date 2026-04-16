
#' Identify divergent loci from tag clusters.
#'
#' @param object A \code{GRanges} object containing the tag clusters.
#' @param ctss A \code{RangedSummarizedExperiment} resulting from
#' \code{CAGEfightR::quantifyCTSSs()}.
#' @param max_gap Maximum distance between tag clusters to be considered
#' part of a divergent locus.
#' @param win_size Size of the flanking windows around the midpoint of the
#' divergent locus.
#' @param inputAssay The assay to use for quantification.
#' 
#' @return A \code{GRanges} object containing the divergent loci.
#' 
#' @export
#' 
#' @importFrom GenomicRanges GRanges findOverlaps flank `seqlevels<-` seqinfo `seqinfo<-` seqlevels strand `strand<-` start end `start<-` `end<-` mcols `mcols<-` names `names<-` subsetByOverlaps isDisjoint swapRanges seqnames
#' @importFrom SummarizedExperiment SummarizedExperiment assay assays assayNames `assayNames<-` rowRanges colData
#' @importFrom IRanges IRanges reduce Views viewSums
#' @importFrom CAGEfightR quantifyClusters calcPooled
#' @importFrom igraph graph_from_edgelist components
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix rowSums
#' @importFrom methods is as
#' 
divergentLoci <- function(object, ctss, max_gap=400, win_size=200, 
                          inputAssay="counts") {
  
  ## Format pooled signal for inputAssay
  object <- CAGEfightR::quantifyClusters(ctss, object)
  object <- suppressWarnings(CAGEfightR::calcPooled(object, inputAssay=inputAssay))
  object <- SummarizedExperiment::rowRanges(object)
  ctss <- suppressWarnings(CAGEfightR::calcPooled(ctss,inputAssay=inputAssay))
  
  
  message("Removing overlapping TCs by strand...")
  ## Split on strand
  TCsByStrand <- CAGEfightR:::splitByStrand(object)
  
  ## Find overlaps between TCs on separate strands
  olaps <- GenomicRanges::findOverlaps(TCsByStrand$'-',TCsByStrand$'+',maxgap=-1,
                        type="any",select="all",ignore.strand=TRUE)
  m_score <-  GenomicRanges::mcols(TCsByStrand$'-')$score
  p_score <-  GenomicRanges::mcols(TCsByStrand$'+')$score
  
  m_rem <- S4Vectors::queryHits(olaps)[which(m_score[S4Vectors::queryHits(olaps)] <= 
                                    p_score[S4Vectors::subjectHits(olaps)])]
  p_rem <- S4Vectors::subjectHits(olaps)[which(p_score[S4Vectors::subjectHits(olaps)] < 
                                      m_score[S4Vectors::queryHits(olaps)])]
  
  ## remove overlapping TCs
  TCsByStrand$'-' <- TCsByStrand$'-'[-m_rem]
  TCsByStrand$'+' <- TCsByStrand$'+'[-p_rem]
  
  message("Finding divergent TC pairs...")
  ## Find divergent TC pairs
  m_pad <- GenomicRanges::flank(TCsByStrand$'-', width=max_gap, start=TRUE, both=FALSE)
  pairs <- GenomicRanges::findOverlaps(m_pad,TCsByStrand$'+',maxgap=-1,type="any",
                        select="all",ignore.strand=TRUE)
  
  ## Find connected components of TC pair graphs
  edge_list <- cbind(names(TCsByStrand$'-')[S4Vectors::queryHits(pairs)],
                     names(TCsByStrand$'+')[S4Vectors::subjectHits(pairs)])
  g <- igraph::graph_from_edgelist(edge_list,directed=FALSE)
  con <- igraph::components(g)
  
  ## Keep only relevant TCs
  object <- object[names(con$membership)]
  
  message("Merging into divergent loci...")
  
  ## Merge connected components into divergent loci
  mergeLoci <- function(minus_pos, plus_pos, minus_score, plus_score)
  {
    m <- length(minus_pos)
    p <- 1
    
    while(minus_pos[m]>plus_pos[p]) {
      
      if (plus_score[p]<minus_score[m])
        p <- p+1
      else
        m <- m-1
    }
    round(mean(c(minus_pos[m],plus_pos[p])))
  }
  
  df <- data.frame(strand=as.character(GenomicRanges::strand(object)),
                   start=GenomicRanges::start(object),
                   end=GenomicRanges::end(object),
                   score=object$score,
                   loci=con$membership,
                   stringsAsFactors=FALSE)
  
  div_mid <- by(df, df$loci, function(d) {
    d <- split(d,d$strand)
    mergeLoci(d$'-'$end,d$'+'$start,d$'-'$score,d$'+'$score)
  })
  
  ## Extract seqnames for loci
  div_chr <- con$membership[match(1:con$no,con$membership)]
  div_chr <- as.character(sapply(names(div_chr), 
                                 function(n) strsplit(n,":")[[1]][1]))
  

  gr <- GenomicRanges::GRanges(seqnames=div_chr,IRanges::IRanges(start=div_mid,end=div_mid))
  GenomicRanges::seqlevels(ctss,pruning.mode="coarse") <- GenomicRanges::seqlevels(gr)
  GenomicRanges::seqinfo(gr) <- GenomicRanges::seqinfo(ctss)

  covByStrand <- CAGEfightR:::splitPooled(methods::as(SummarizedExperiment::rowRanges(ctss),
                                          "GRanges"))
  
  cat("\r")
  message("Calculating directionality...")
  
  win_1 <- GenomicRanges::flank(gr,width=win_size,start=TRUE,both=FALSE)
  win_2 <- GenomicRanges::flank(gr,width=win_size,start=FALSE,both=FALSE)
  
  ## Quantify strand-wise in flanking windows around midpoint
  M1 <- unlist(IRanges::viewSums(IRanges::Views(covByStrand$`-`, win_1)))
  P2 <- unlist(IRanges::viewSums(IRanges::Views(covByStrand$`+`, win_2)))
  M2 <- unlist(IRanges::viewSums(IRanges::Views(covByStrand$`-`, win_2)))
  P1 <- unlist(IRanges::viewSums(IRanges::Views(covByStrand$`+`, win_1)))
  
  ## Calculate directionality
  pooled_directionality <- (P2-M1) / (P2+M1)
  
  ## Test if divergent
  divergent <- (M1>P1) & (P2>M2)
  
  message("Calculating coverage across samples...")
  
  ## Quantify strand-wise in flanking windows around midpoint
  GenomicRanges::strand(win_1) <- "-"
  GenomicRanges::strand(win_2) <- "+"
  mat_2_plus <- 
    suppressMessages(SummarizedExperiment::assay(CAGEfightR::quantifyClusters(
      ctss, win_2, 
      inputAssay = inputAssay),
      inputAssay) > 0)
  mat_1_minus <- 
    suppressMessages(SummarizedExperiment::assay(CAGEfightR::quantifyClusters(
      ctss, win_1, 
      inputAssay = inputAssay),
      inputAssay) > 0)
  
  ## Quntify number of bidirectional cases (both strands expressed)
  bidirectional <- Matrix::rowSums(mat_1_minus & mat_2_plus)
  
  message("Preparing output...")
  
  ## Build GRanges object
  start(gr) <- GenomicRanges::start(gr)-win_size
  end(gr) <- GenomicRanges::end(gr)+win_size
  gr$score <- M1+P2
  gr$thick <- IRanges::IRanges(start=div_mid,width=1)
  
  GenomicRanges::mcols(gr)[, "directionality"] <- pooled_directionality
  GenomicRanges::mcols(gr)[, "bidirectionality"] <- bidirectional
  
  ids <- paste0(GenomicRanges::seqnames(gr), ":", GenomicRanges::start(gr), "-", GenomicRanges::end(gr))
  names(gr) <- ids
  names(gr$thick) <- ids
  
  ## Remove non-divergent cases
  gr[divergent]
}

#' Quantify the expression of divergent loci (identified by 
#' \code{divergentLoci}).
#'
#' @param loci A \code{GRanges} object containing the divergent loci.
#' @param ctss A \code{RangedSummarizedExperiment} resulting from
#' \code{CAGEfightR::quantifyCTSSs()}.
#' @param inputAssay The assay to use for quantification.
#' @param requireDisjoint Logical indicating whether to require divergent 
#' loci to be disjoint.
#' 
#' @return A \code{SummarizedExperiment} object containing the quantified
#' divergent loci.
#' 
#' @export
#' 
#' @importFrom S4Vectors SimpleList
quantifyDivergentLoci <- function(loci, ctss, inputAssay="counts", 
                                  requireDisjoint=TRUE) {
  
  res <- quantifyStrandwiseDivergentLoci(loci, ctss, 
                                         inputAssay, 
                                         requireDisjoint)
  
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

#' Quantify the expression of disjoint tag clusters.
#' 
#' @param object A \code{RangedSummarizedExperiment} resulting from
#' \code{CAGEfightR::quantifyCTSSs()}.
#' @param clusters A \code{GRanges} object containing the tag clusters.
#' @param inputAssay The assay to use for quantification.
#' 
#' @return A \code{SummarizedExperiment} object containing the quantified
#' tag clusters.
#' 
#' @export
#' 
#' @importFrom BiocParallel bplapply
quantifyClustersOlap <- function(object, clusters, inputAssay) {
  
  revmap <- IRanges::reduce(clusters, with.revmap=TRUE)$revmap
  max.nr <- max(sapply(revmap, length))
  ids <- lapply(1:max.nr, function(level) {
    setdiff(sapply(revmap, `[`, level),NA)})
  
  res <- BiocParallel::bplapply(1:max.nr, function(i) {
    clu <- clusters[ids[[i]]]
    obj <- object
    GenomicRanges::seqinfo(obj) <- GenomicRanges::seqinfo(clu)
    
    ## special case, CAGEfightR::quantifyClusters does not deal with 
    ## granges objects of length 1
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
  
  o <- SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(mat),
    rowRanges = clusters,
    colData = SummarizedExperiment::colData(object))
  
  SummarizedExperiment::assayNames(o) <- inputAssay
  
  o
}

## Helper function, quantify the expression of divergent loci in a 
## strand-wise manner

#' import BiocGenerics
#' 
#' @param loci A \code{GRanges} object containing the divergent loci.
#' @param ctss A \code{RangedSummarizedExperiment} resulting from
#' \code{CAGEfightR::quantifyCTSSs()}.
#' @param inputAssay The assay to use for quantification.
#' @param requireDisjoint Logical indicating whether to require divergent
#' loci to be disjoint.
quantifyStrandwiseDivergentLoci <- function(loci, ctss, inputAssay = "counts", 
                                            requireDisjoint=TRUE) {
  win_1 <- loci
  end(win_1) <- start(loci$thick) - 1
  strand(win_1) <- "-"
  
  win_2 <- loci
  start(win_2) <- start(loci$thick) + 1
  strand(win_2) <- "+"
  
  m1 <- matrix()
  if (!requireDisjoint && !GenomicRanges::isDisjoint(win_1)) {
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
  if (!requireDisjoint && !GenomicRanges::isDisjoint(win_2)) {
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
  
  res <- list(m,p)
  names(res) <- c("-","+")
  
  res
}



#' Identify divergent loci from tag clusters with summit centering.
#'
#' @param object A \code{GRanges} object containing the tag clusters.
#' @param ctss A \code{RangedSummarizedExperiment} resulting from
#' \code{CAGEfightR::quantifyCTSSs()}.
#' @param max_gap Maximum distance between tag clusters to be considered
#' part of a divergent locus.
#' @param win_size Size of the flanking windows around the midpoint of the
#' divergent locus.
#' @param inputAssay The assay to use for quantification.
#' 
#' @return A \code{GRanges} object containing the divergent loci.
#' 
#' @export
#' 
#' @importFrom GenomicRanges GRanges findOverlaps flank seqinfo `seqinfo<-` seqlevels `seqlevels<-` start end `start<-` `end<-` strand `strand<-` mcols `mcols<-` names `names<-` seqnames subsetByOverlaps swapRanges
#' @importFrom SummarizedExperiment assay assays assayNames `assayNames<-` rowRanges colData
#' @importFrom IRanges IRanges reduce Views viewSums
#' @importFrom CAGEfightR quantifyClusters calcPooled
#' @importFrom data.table data.table setorder
#' @importFrom igraph graph_from_edgelist components
#' @importFrom BiocParallel bplapply
#' 
divergentLociSummit<- function(object, ctss, max_gap=400, win_size=200, 
                                  inputAssay="counts") {

  ## Format pooled signal for inputAssay
  object <- CAGEfightR::quantifyClusters(ctss, object)
  object <- suppressWarnings(CAGEfightR::calcPooled(object, inputAssay=inputAssay))
  object <- SummarizedExperiment::rowRanges(object)
  ctss <- suppressWarnings(CAGEfightR::calcPooled(ctss,inputAssay=inputAssay))
  
  message("Removing overlapping TCs by strand...")
  object <- GenomicRanges::swapRanges(object) #Summit focused
  TCsByStrand <- CAGEfightR:::splitByStrand(object)
  
  ## Find overlapping and book-ended summits between strands
  olaps <- GenomicRanges::findOverlaps(TCsByStrand$'-',TCsByStrand$'+',maxgap=0,type="any",select="all",ignore.strand=TRUE)
  m_score <-  GenomicRanges::mcols(TCsByStrand$'-')$score
  p_score <-  GenomicRanges::mcols(TCsByStrand$'+')$score
  
  m_rem <- S4Vectors::queryHits(olaps)[which(m_score[S4Vectors::queryHits(olaps)] <= p_score[S4Vectors::subjectHits(olaps)])]
  p_rem <- S4Vectors::subjectHits(olaps)[which(p_score[S4Vectors::subjectHits(olaps)] < m_score[S4Vectors::queryHits(olaps)])]
  
  ## remove overlapping TCs
  if (length(m_rem)>0) {
    TCsByStrand$'-' <- TCsByStrand$'-'[-m_rem]
  }
  if (length(p_rem)>0) {
    TCsByStrand$'+' <- TCsByStrand$'+'[-p_rem]
  }
  
  message("Finding divergent TC pairs...")
  ## Find divergent TC pairs
  m_pad <- GenomicRanges::flank(TCsByStrand$'-', width=max_gap, start=TRUE, both=FALSE)
  pairs <- GenomicRanges::findOverlaps(m_pad,TCsByStrand$'+',maxgap=-1,type="any",select="all",ignore.strand=TRUE)
  
  ## Find connected components of TC pair graphs
  edge_list <- cbind(names(TCsByStrand$'-')[S4Vectors::queryHits(pairs)],
                     names(TCsByStrand$'+')[S4Vectors::subjectHits(pairs)])
  g <- igraph::graph_from_edgelist(edge_list,directed=FALSE)
  con <- igraph::components(g)
  
  ## Keep only relevant TCs
  object <- object[names(con$membership)]
  
  message("Merging into divergent loci...")
  covByStrand <- CAGEfightR:::splitPooled(methods::as(SummarizedExperiment::rowRanges(ctss),"GRanges"))
  
  mergeLoci <- function(m, p) {
    
    m.i <- 1
    p.i <- 1
    
    while ((p[p.i,start] < m[m.i,start]) || ((p[p.i,start]-m[m.i,start]) > max_gap)) {
      if ((p[p.i,score] < m[m.i,score]) && (p.i < nrow(p)) && ((p[p.i+1,start]-m[m.i,start]) < max_gap)) {
        p.i <- p.i + 1
      } else {
        m.i <- m.i + 1
      }
    }
    floor(m[m.i,start] + ((p[p.i,start] - m[m.i,start]) / 2))
  }
  
  dt <- data.table::data.table(strand=as.character(GenomicRanges::strand(object)),
                   start=GenomicRanges::start(object),
                   score=object$score,
                   locus=con$membership,
                   stringsAsFactors=FALSE)
  
  dt <- split(dt, dt$strand)
  data.table::setorder(dt$'-', -score, start) ## order by score, then by position (prioritize rightmost)
  data.table::setorder(dt$'+', -score, -start) ## order by score, then by position (prioritize leftmost)
  
  loci.m <- split(dt$'-', dt$'-'$'locus')
  loci.p <- split(dt$'+', dt$'+'$'locus')
  
  mid <- BiocParallel::bplapply(1:length(loci.m), function(i) mergeLoci(loci.m[[i]], loci.p[[i]]))
  mid <- unlist(mid)
  
  ## Extract seqnames for loci
  div_chr <- con$membership[match(1:con$no,con$membership)]
  div_chr <- as.character(sapply(names(div_chr), function(n) strsplit(n,":")[[1]][1]))
  
  
  gr <- GenomicRanges::GRanges(seqnames=div_chr,IRanges::IRanges(start=mid,end=mid))
  GenomicRanges::seqlevels(ctss,pruning.mode="coarse") <- GenomicRanges::seqlevels(gr)
  GenomicRanges::seqinfo(gr) <- GenomicRanges::seqinfo(ctss)
  
  cat("\r")
  message("Calculating directionality...")
  
  win_1 <- GenomicRanges::flank(gr,width=win_size,start=TRUE,both=FALSE)
  win_2 <- GenomicRanges::flank(gr,width=win_size,start=FALSE,both=FALSE)
  
  covByStrand <- CAGEfightR:::splitPooled(methods::as(SummarizedExperiment::rowRanges(ctss),"GRanges"))
  
  ## Quantify strand-wise in flanking windows around midpoint
  M1 <- unlist(IRanges::viewSums(IRanges::Views(covByStrand$`-`, win_1)))
  P2 <- unlist(IRanges::viewSums(IRanges::Views(covByStrand$`+`, win_2)))
  M2 <- unlist(IRanges::viewSums(IRanges::Views(covByStrand$`-`, win_2)))
  P1 <- unlist(IRanges::viewSums(IRanges::Views(covByStrand$`+`, win_1)))
  
  ## Calculate directionality
  pooled_directionality <- (P2-M1) / (P2+M1)
  
  ## Test if divergent
  divergent <- (M1>P1) & (P2>M2)
  
  message("Calculating coverage across samples...")
  
  ## Quantify strand-wise in flanking windows around midpoint
  GenomicRanges::strand(win_1) <- "-"
  GenomicRanges::strand(win_2) <- "+"
  mat_2_plus <- suppressMessages(SummarizedExperiment::assay(CAGEfightR::quantifyClusters(ctss, win_2, inputAssay = inputAssay),inputAssay) > 0)
  mat_1_minus <- suppressMessages(SummarizedExperiment::assay(CAGEfightR::quantifyClusters(ctss, win_1, inputAssay = inputAssay),inputAssay) > 0)
  
  ## Quntify number of bidirectional cases (both strands expressed)
  bidirectional <- Matrix::rowSums(mat_1_minus & mat_2_plus)
  
  message("Preparing output...")
  
  ## Build GRanges object
  start(gr) <- GenomicRanges::start(gr)-win_size
  end(gr) <- GenomicRanges::end(gr)+win_size
  gr$score <- M1+P2
  gr$thick <- IRanges::IRanges(start=mid,width=1)
  
  GenomicRanges::mcols(gr)[, "directionality"] <- pooled_directionality
  GenomicRanges::mcols(gr)[, "bidirectionality"] <- bidirectional
  GenomicRanges::mcols(gr)[, "divergent"] <- divergent
  
  ids <- paste0(GenomicRanges::seqnames(gr), ":", GenomicRanges::start(gr), "-", GenomicRanges::end(gr))
  names(gr) <- ids
  names(gr$thick) <- ids
  
  ## Remove non-divergent cases
  gr
}

#'
#' Identify divergent loci from tag clusters with summit centering. The function is similar to \code{divergentLociSummit} 
#' but does not require tag clusters to be passed to the function. This is handeled internally by calling 
#' \code{CAGEfightR::clusterUnidirectionally} on the pooled CTSS signal.
#' 
#' @param ctss A \code{RangedSummarizedExperiment} resulting from
#' \code{CAGEfightR::quantifyCTSSs()}.
#' @param max_gap Maximum distance between tag clusters to be considered
#' part of a divergent locus.
#' @param win_size Size of the flanking windows around the midpoint of the
#' divergent locus.
#' @param inputAssay The assay to use for quantification.
#' @param callingAssay The assay to use for calling TCs.
#' 
#' @return A \code{GRanges} object containing the divergent loci.
#' 
#' @export
#' 
#' @importFrom data.table data.table setorder
#' @importFrom igraph graph_from_edgelist components
#' @importFrom BiocParallel bplapply
#' 
divergentLociTCsSummit <- function(ctss, max_gap = 400, win_size = 200, inputAssay = "counts", callingAssay = "counts") {
    
    base::message(base::paste("Running DL calling for",base::colnames(ctss)))
    
    ## Pool & subset
    ctss.o <- ctss
    ctss <- base::suppressWarnings(
        CAGEfightR::calcPooled(ctss, inputAssay = callingAssay)
      )
    ctss <- base::subset(ctss, score > 0)
    
    #Call TCs for each sample
    object <- CAGEfightR::clusterUnidirectionally(ctss)

    #Call divergent loci with summit centering
    gr <- divergentLociSummit(object, ctss.o, max_gap = max_gap, win_size = win_size, inputAssay = inputAssay)

    if(base::length(gr) > 0) {
        gr$sample <- base::colnames(ctss)
    }
    
    return(gr)
  }