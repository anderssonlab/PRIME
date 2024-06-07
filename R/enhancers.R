require("CAGEfightR")
source("CAGEfightR_extensions/utils.R")

## Object: GRanges
## ctss: RangedSummarizedExperiment
divergentLoci <- function(object, ctss, max_gap=400, win_size=200, inputAssay="counts") {

    message("Removing overlapping TCs by strand...")
    ## Split on strand
    TCsByStrand <- splitByStrand(object)

    ## Find overlaps between TCs on separate strands
    olaps <- findOverlaps(TCsByStrand$'-',TCsByStrand$'+',maxgap=-1,type="any",select="all",ignore.strand=TRUE)
    m_score <-  mcols(TCsByStrand$'-')$score
    p_score <-  mcols(TCsByStrand$'+')$score
    
    m_rem <- queryHits(olaps)[which(m_score[queryHits(olaps)] <= p_score[subjectHits(olaps)])]
    p_rem <- subjectHits(olaps)[which(p_score[subjectHits(olaps)] < m_score[queryHits(olaps)])]

    ## remove overlapping TCs
    TCsByStrand$'-' <- TCsByStrand$'-'[-m_rem]
    TCsByStrand$'+' <- TCsByStrand$'+'[-p_rem]

    message("Finding divergent TC pairs...")
    ## Find divergent TC pairs
    m_pad <- flank(TCsByStrand$'-', width=max_gap, start=TRUE, both=FALSE)
    pairs <- findOverlaps(m_pad,TCsByStrand$'+',maxgap=-1,type="any",select="all",ignore.strand=TRUE)

    ## Find connected components of TC pair graphs
    edge_list <- cbind(names(TCsByStrand$'-')[queryHits(pairs)],
                     names(TCsByStrand$'+')[subjectHits(pairs)])
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

    df <- data.frame(strand=as.character(strand(object)),
                     start=start(object),
                     end=end(object),
                     score=object$score,
                     loci=con$membership,
                     stringsAsFactors=FALSE)
    
    div_mid <- by(df, df$loci, function(d) {
        d <- split(d,d$strand)
        mergeLoci(d$'-'$end,d$'+'$start,d$'-'$score,d$'+'$score)
    })
            
    ## Extract seqnames for loci
    div_chr <- con$membership[match(1:con$no,con$membership)]
    div_chr <- as.character(sapply(names(div_chr), function(n) strsplit(n,":")[[1]][1]))
    
    covByStrand <- splitPooled(methods::as(rowRanges(ctss),"GRanges"))
    gr <- GRanges(seqnames=div_chr,IRanges(start=div_mid,end=div_mid))
    seqinfo(gr) <- seqinfo(ctss)

    cat("\r")
    message("Calculating directionality...")

    win_1 <- flank(gr,width=win_size,start=TRUE,both=FALSE)
    win_2 <- flank(gr,width=win_size,start=FALSE,both=FALSE)

    ## Quantify strand-wise in flanking windows around midpoint
    M1 <- unlist(viewSums(Views(covByStrand$`-`, win_1)))
    P2 <- unlist(viewSums(Views(covByStrand$`+`, win_2)))
    M2 <- unlist(viewSums(Views(covByStrand$`-`, win_2)))
    P1 <- unlist(viewSums(Views(covByStrand$`+`, win_1)))

    ## Calculate directionality
    pooled_directionality <- (P2-M1) / (P2+M1)

    ## Test if divergent
    divergent <- (M1>P1) & (P2>M2)
    
    message("Calculating coverage across samples...")

    ## Quantify strand-wise in flanking windows around midpoint
    strand(win_1) <- "-"
    strand(win_2) <- "+"
    mat_2_plus <- suppressMessages(assay(quantifyClusters(ctss, win_2, inputAssay = inputAssay),inputAssay) > 0)
    mat_1_minus <- suppressMessages(assay(quantifyClusters(ctss, win_1, inputAssay = inputAssay),inputAssay) > 0)

    ## Quntify number of bidirectional cases (both strands expressed)
    bidirectional <- rowSums(mat_1_minus & mat_2_plus)

    message("Preparing output...")

    ## Build GRanges object
    start(gr) <- start(gr)-win_size
    end(gr) <- end(gr)+win_size
    gr$score <- M1+P2
    gr$thick <- IRanges(start=div_mid,width=1)
    
    mcols(gr)[, "directionality"] <- pooled_directionality
    mcols(gr)[, "bidirectionality"] <- bidirectional

    ids <- paste0(seqnames(gr), ":", start(gr), "-", end(gr))
    names(gr) <- ids
    names(gr$thick) <- ids
    
    ## Remove non-divergent cases
    gr[divergent]
}
