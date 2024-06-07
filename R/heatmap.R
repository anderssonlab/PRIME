require("CAGEfightR")
require("assertthat")

source("CAGEfightR_extensions/utils.R")

## data: GRanges with value to consider in specified column (score)
## regions: GRanges with regions to extract data for (must be of the same size)
## transform: function applied on data for each region (e.g. to reduce dimension=, defaults to identity

heatmapData <- function(regions, data, column="score", transform_fn=identity, ...) {
    
    ## Check regions
    assert_that(length(unique(width(regions)))==1,
                column %in% colnames(mcols(data)))

    sl <- intersect(seqlevels(regions),seqlevels(data))
    if (!all(seqlevels(regions) %in% sl) || !all(seqlevels(data) %in% sl)) {
        warning("seqlevels differ between regions and data GRanges obejects, subsetting to intersection")
        seqlevels(regions, pruning.mode="coarse") <- seqlevels(data, pruning.mode="coarse") <- sl
    }   
    
    regionsByStrand <- splitByStrand(regions)
    dataByStrand <- splitPooled(data, weight=column)

    nr <- names(regionsByStrand)[sapply(regionsByStrand, function(x) length(x)>0)]
    nd <- names(dataByStrand)[sapply(dataByStrand, function(x) any(sapply(x,function(y) any(y!=0))))]
    
    res <- lapply(nr, function(r) {
        message("extracting data for strand: ", r)
        
        dat <- lapply(nd, function(d) {
            message("   ",d)
            vl <- Views(dataByStrand[[d]], regionsByStrand[[r]])
            n <- paste0(unlist(lapply(names(vl), function(nv) rep(nv,length(vl[[nv]])))), ":",
                        unlist(lapply(vl,start)), "-", unlist(lapply(vl,end)), ";", d)
            
            vlen <- sapply(vl,length)
            
            m <- do.call("rbind",lapply(vl[vlen>0], function(v) t(viewApply(v,function(x) transform_fn(as.vector(x,mode="numeric"),...)))))
            rownames(m) <- n
            m
        })
        names(dat) <- nd

        dat
    })

    names(res) <- nr

    res
}
