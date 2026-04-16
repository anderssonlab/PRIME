#' Export bigWig files for each replicate, optionally merging strands.
#'
#' @param object A \code{RangedSummarizedExperiment} object containing the CTSSs.
#' @param dir A character string of the directory to export the bigWig files to.
#' @param replicates A character vector of replicate names to export. Default is "all".
#' @param inputAssay A character string of the assay to use for the bigWig scores. Default is "TPM".
#' @param splitByStrand A logical indicating whether to export strands separately. Default is TRUE.
#'
#' @return No return value. The function exports bigWig files to the specified directory.
#'
#' @importFrom assertthat assert_that
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom rtracklayer export
#' @importFrom SummarizedExperiment assay rowRanges colData
#' @importFrom GenomicRanges GRanges
#' @importFrom BiocGenerics strand start end
#' @importFrom IRanges reduce
#' @importFrom S4Vectors mcols aggregate
#' @importFrom magrittr %>%
#' @export
#'
writeBw <- function(object, dir, replicates = "all", inputAssay = "TPM", splitByStrand = TRUE) {
  if (replicates == "all") {
    replicates <- base::rownames(SummarizedExperiment::colData(object))
  }

  ## Check consistency of replicates
  assertthat::assert_that(
    base::length(SummarizedExperiment::rowRanges(object)) == base::nrow(SummarizedExperiment::assay(object, inputAssay)),
    base::all(replicates %in% base::rownames(SummarizedExperiment::colData(object)))
  )
  base::stopifnot("Chosen directory doesn't exist!" = base::dir.exists(dir))

  gr <- GenomicRanges::GRanges(SummarizedExperiment::rowRanges(object))

  ## Parallel export across all chosen replicates
  foreach::foreach(i = 1:base::length(replicates)) %dopar% {
    ## Print replicate when exiting
    base::on.exit(base::cat(base::paste0("Exit due to ", replicates[i], ".\n")), add = TRUE)

    ## Convert individual replicate to GRanges object
    S4Vectors::mcols(gr)["score"] <- SummarizedExperiment::assay(object, inputAssay)[, replicates[i]]

    ## Export replicate's plus and minus strands separately
    if (splitByStrand) {
      for (s in base::c("+", "-")) {
        grs <- gr[BiocGenerics::strand(gr) == s, ]
        s <- base::ifelse(s == "+", "plus", "minus")

        fn <- base::paste0(dir, replicates[i], ".", inputAssay, ".", s, ".bw")
        rtracklayer::export(
          object = grs,
          con = fn,
          format = "bigWig"
        )
        base::message(base::paste(fn, "done.", sep = " "))
      }
    } else {
      # Aggregate counts from both strands
      grred <- IRanges::reduce(gr, ignore.strand = TRUE, with.revmap = TRUE, min.gapwidth = 0)
      S4Vectors::mcols(grred) <- S4Vectors::aggregate(
        gr, S4Vectors::mcols(grred)$revmap,
        score = base::sum(score), drop = FALSE
      )

      # Keep only score metadata column
      grred@elementMetadata <- grred@elementMetadata[, "score", drop = FALSE]

      fn <- base::paste0(dir, replicates[i], ".", inputAssay, ".bw")
      rtracklayer::export(
        object = grred,
        con = fn,
        format = "bigWig"
      )
      base::message(base::paste(fn, "done.", sep = " "))
    }
  }
}

#' Read bigWig files for a specified genomic range and return a tibble of scores.
#'
#' @param files A character vector of bigWig file paths to read.
#' @param range A \code{GRanges} object specifying the genomic range to read from the bigWig files.
#'
#' @return A tibble of scores for the specified genomic range.
#'
#' @importFrom rtracklayer BigWigSelection
#' @importFrom rtracklayer import
#' @importFrom purrr map
#' @importFrom dplyr bind_rows mutate case_when
#' @importFrom stringr str_replace
#' @importFrom assertthat assert_that
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom methods is
#'
#' @export
#'
readRange <- function(files, range) {
  assertthat::assert_that(
    is.character(files),
    length(files) >= 1,
    all(file.exists(files)),
    methods::is(range, "GRanges"),
    length(range) >= 1
  )

  ## Prepare range for rtracklayer
  ran_sel <- rtracklayer::BigWigSelection(
    ranges = range
  )
  ## Reading in scores in range from files serially
  raw <- purrr::map(
    .x = files,
    .f = function(f) {
      s <- base::basename(f) %>%
        stringr::str_replace(., "^(.*?).bw$", "\\1")
      return(
        base::suppressWarnings(
          rtracklayer::import(
            con = f,
            as = "GRanges",
            selection = ran_sel
          ) %>%
            tibble::as_tibble() %>%
            dplyr::mutate(
              sample = s,
              strand = stringr::str_replace(
                s, "^(.*)\\.(.*?)$", "\\2"
              ),
              score = dplyr::case_when(
                strand == "minus" ~ (-1) * score,
                TRUE ~ score
              )
            )
        )
      )
    }
  ) %>%
    dplyr::bind_rows()

  ## Fill missing scores across samples with zeros
  out <- fillGaps(raw)

  return(out)
}

#' Fill in missing scores across samples with zeros. Internal function used by \code{readRange}.
#'
#' @param data A tibble of scores generated by \code{readRange}.
#'
#' @return A tibble with missing scores filled in with zeros.
#'
#' @importFrom dplyr distinct pull across
#' @importFrom tidyr pivot_wider replace_na pivot_longer
#' @importFrom stringr str_replace
#' @importFrom tidyselect all_of
#'
fillGaps <- function(data) {
  assertthat::assert_that(
    is.data.frame(data),
    all(c("sample", "score") %in% colnames(data))
  )

  columns <- data %>%
    dplyr::distinct(sample) %>%
    dplyr::pull(sample)

  df_filled <- data %>%
    tidyr::pivot_wider(
      names_from = "sample",
      values_from = "score"
    ) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = columns,
        .fns = function(x) {
          tidyr::replace_na(x, 0)
        }
      )
    ) %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(columns),
      names_to = "sample",
      values_to = "score"
    ) %>%
    dplyr::mutate(
      sample = stringr::str_replace(
        sample, ".minus|.plus", ""
      )
    )
  return(df_filled)
}

#' Plot tracks for a specified genomic range using ggplot2.
#'
#' @param data A tibble of scores for the specified genomic range, generated by \code{readRange}.
#' @param range A \code{GRanges} object specifying the genomic range to plot
#'
#' @return A ggplot object representing the tracks for the specified genomic range.
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_x_continuous scale_color_manual scale_linetype_manual facet_wrap theme_bw theme element_blank element_text
#' @importFrom grid unit
#' @importFrom BiocGenerics start end
#' @importFrom methods is
#'
#' @export
#'
plotTracks <- function(data, range) {
  assertthat::assert_that(
    is.data.frame(data),
    all(c("start", "score", "strand", "sample") %in% colnames(data)),
    methods::is(range, "GRanges"),
    length(range) >= 1
  )

  start_pos <- BiocGenerics::start(range)
  end_pos <- BiocGenerics::end(range)

  plot <- ggplot2::ggplot(
    data = data,
    mapping = ggplot2::aes(
      x = start,
      y = score,
      colour = strand,
      linetype = strand
    )
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      size = 0.5
    ) +
    ggplot2::scale_x_continuous(
      limits = base::c(start_pos, end_pos),
      position = "top"
    ) +
    ggplot2::scale_color_manual(
      limits = base::c("minus", "plus"),
      values = base::c("#FF6400", "#632770")
    ) +
    ggplot2::scale_linetype_manual(
      limits = base::c("minus", "plus"),
      breaks = base::c("minus", "plus"),
      values = base::c("solid", "solid")
    ) +
    ggplot2::facet_wrap(
      facets = . ~ sample,
      ncol = 1,
      strip.position = "right"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.05, "lines"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = grid::unit(base::c(0, 0, 0, 0), "cm")
    )
  return(plot)
}