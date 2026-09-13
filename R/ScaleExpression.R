#' Centre and scale features within a Bioconductor experiment
#'
#' Each feature is centred and/or scaled across pixels using base R's
#' `scale()`. Zero-variance features become zero after centring and scaling.
#' The input assay is preserved and the result is stored in `scaled` by
#' default. This materialises a dense matrix; use implicit centring in PCA
#' instead when a full scaled matrix would exceed available memory.
#'
#' @param data A SingleCellExperiment, including SpatialExperiment.
#' @param assay Primary (`"main"`) or alternative experiment name.
#' @param slot Input expression assay, usually `"logcounts"`.
#' @param center,scale Logical values passed to base `scale()`.
#' @param outputAssay Name of a new assay for scaled values. Must differ from
#'   the input assay and must not already exist.
#'
#' @return The input container with a scaled assay in the selected experiment.
#' @export
#' @examples
#' counts <- matrix(seq_len(20), nrow = 4,
#'     dimnames = list(paste0("mz-", 101:104), paste0("p", 1:5)))
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
#' sce <- normalizeSMData(sce, "LogNormalize", verbose = FALSE)
#' sce <- scaleSMData(sce)
#' SummarizedExperiment::assayNames(sce)
scaleSMData <- function(data, assay = "main", slot = "logcounts",
                        center = TRUE, scale = TRUE, outputAssay = "scaled") {
  .requireExperiment(data, "SingleCellExperiment")
  for (flag in list(center, scale)) {
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) {
      stop("center and scale must each be TRUE or FALSE.", call. = FALSE)
    }
  }
  experiment <- .experimentForAssay(data, assay)
  if (length(outputAssay) != 1L || is.na(outputAssay) || !nzchar(outputAssay) ||
      outputAssay %in% SummarizedExperiment::assayNames(experiment)) {
    stop("outputAssay must name a new assay; existing values are not overwritten.",
         call. = FALSE)
  }
  values <- as.matrix(.assayData(data, assay, slot))
  if (!nrow(values) || ncol(values) < 2L || any(!is.finite(values))) {
    stop("Scaling requires finite values and at least two pixels.", call. = FALSE)
  }
  scaled <- base::scale(t(values), center = center, scale = scale)
  centers <- attr(scaled, "scaled:center")
  scales <- attr(scaled, "scaled:scale")
  if (isTRUE(scale)) {
    constant <- which(scales == 0)
    scaled[, constant] <- 0
  }
  if (any(!is.finite(scaled))) {
    stop("Scaling produced non-finite values; check the intensity range.", call. = FALSE)
  }
  SummarizedExperiment::assay(experiment, outputAssay) <- t(scaled)
  S4Vectors::metadata(experiment)$spamtp_scaling <- list(
    input_assay = slot, output_assay = outputAssay,
    center = centers, scale = scales)
  .replaceExperiment(data, experiment, assay)
}
