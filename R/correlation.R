source("CAGEfightR_extensions/heatmap.R")

## Calculates the auto correlation of CTSSs over tag clusters according to pooled values
## and the cross correlation with opposite strand CTSSs

## object: GRanges
## pooled: RangedSummarizedExperiment
corrprofiles <- function(object, pooled, flank=300, lag.max=300) {

    ## flank TCs around their summits
    start(object) <- start(object$thick)-flank
    end(object) <- start(object$thick)+flank

    ## retrieve data around regions
    message("Retrieving data")
    data <- heatmapData(object, as(rowRanges(pooled),"GRanges"))

    strand <- as.character(strand(object))
    sense_matrix <- matrix(0, ncol=2*flank+1, nrow=length(object))
    antisense_matrix <- matrix(0, ncol=2*flank+1, nrow=length(object))

    sense_matrix[strand=="+",] <- data[["+"]][["+"]]
    sense_matrix[strand=="-",] <- data[["-"]][["-"]][,(2*flank+1):1]
    antisense_matrix[strand=="+",] <- data[["+"]][["-"]]
    antisense_matrix[strand=="-",] <- data[["-"]][["+"]][,(2*flank+1):1]

    message("Calculating autocorrelation profiles for sense strand")
    ac_matrix <- t(sapply(1:length(object), function(i) acf(sense_matrix[i,],plot=FALSE,lag.max=lag.max)$acf[,1,1]))

    message("Calculating crosscorrelation profiles")
    cc_matrix <- t(sapply(1:length(object), function(i) ccf(sense_matrix[i,],antisense_matrix[i,],plot=FALSE,lag.max=lag.max)$acf[,1,1]))

    rownames(ac_matrix) <- names(object)
    rownames(cc_matrix) <- names(object)

    return(list(ac=ac_matrix,cc=cc_matrix))
}
