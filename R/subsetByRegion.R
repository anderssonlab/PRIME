require("CAGEfightR")

## gr: GRanges
## ctss: RangedSummarizedExperiment
subsetByRegion <- function(gr, ctss) {

    ## Find overlaps
    message("Finding overlaps...")

    hits <- findOverlaps(query = ctss,
                         subject = gr,
                         select = "arbitrary")
    hits <- factor(hits, levels=seq_along(gr))
    ids <- as.numeric(unlist(split(1:length(hits), hits)))

    ctss[ids,]
}
