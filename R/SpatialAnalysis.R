#' Find features correlated with a metabolite, gene, or pixel group
#'
#' Computes Pearson correlations across matched pixels. Genes are read from
#' altExp without assigning artificial m/z values. Constant features return NA.
#' @param data A SpatialExperiment or aligned MSImagingExperiment.
#' @param mz Numeric m/z or exact primary feature name.
#' @param gene Exact feature name in ST.assay.
#' @param ident A colData column defining groups. Supply exactly one of mz,
#'   gene and ident. Groups are compared with binary membership vectors.
#' @param SM.assay Primary assay or modality; `main` selects the primary
#'   experiment. Legacy `SPM` and `Spatial` also select the primary experiment.
#' @param ST.assay Alternative experiment name, or NULL for metabolites only.
#' @param SM.slot Assay containing metabolite values.
#' @param ST.slot Assay containing transcriptomic values.
#' @param nfeatures Maximum results per reference; NULL returns all features.
#' @return A data frame with features, correlation, modality, mz, rank and,
#'   for multiple reference groups, ident.
#' @export
#' @examples
#' x <- SpatialExperiment::SpatialExperiment(
#'   assays = list(counts = rbind(a = 1:4, b = 4:1)),
#'   rowData = S4Vectors::DataFrame(mz = c(100, 200)),
#'   spatialCoords = cbind(x = 1:4, y = 0),
#'   colData = S4Vectors::DataFrame(row.names = paste0("p", 1:4)))
#' findCorrelatedFeatures(x, mz = 100)
findCorrelatedFeatures <- function(
    data, mz = NULL, gene = NULL, ident = NULL, SM.assay = "main",
    ST.assay = NULL, SM.slot = "counts", ST.slot = "counts", nfeatures = 10
) {
  data <- .nativeSpatialObject(data)
  if (sum(!vapply(list(mz, gene, ident), is.null, logical(1))) != 1L)
    stop("Supply exactly one of mz, gene, and ident.", call. = FALSE)
  .validateFeatureCount(nfeatures, "nfeatures", allowNull = TRUE)
  metabolites <- .assayData(data, SM.assay, SM.slot)
  matrices <- list(metabolite = metabolites)
  if (!is.null(ST.assay)) {
    if (!ST.assay %in% SingleCellExperiment::altExpNames(data))
      stop("ST.assay must name an alternative experiment.", call. = FALSE)
    transcripts <- .assayData(data, ST.assay, ST.slot)
    if (is.null(colnames(transcripts)) ||
        !setequal(colnames(transcripts), colnames(metabolites)))
      stop("Modalities must contain the same named pixels.", call. = FALSE)
    matrices$gene <- transcripts[, colnames(metabolites), drop = FALSE]
  }
  if (!is.null(mz)) {
    feature <- if (is.numeric(mz)) findNearestMZ(data, mz, SM.assay) else mz
    if (length(feature) != 1L || !feature %in% rownames(metabolites))
      stop("mz must identify one primary feature.", call. = FALSE)
    references <- list(reference = as.numeric(metabolites[feature, ]))
  } else if (!is.null(gene)) {
    if (is.null(ST.assay) || length(gene) != 1L ||
        !gene %in% rownames(matrices$gene))
      stop("gene must identify one feature in ST.assay.", call. = FALSE)
    references <- list(reference = as.numeric(matrices$gene[gene, ]))
  } else {
    metadata <- SummarizedExperiment::colData(data)
    if (length(ident) != 1L || !ident %in% colnames(metadata))
      stop("ident must name a colData column.", call. = FALSE)
    groups <- as.character(metadata[[ident]])
    if (anyNA(groups)) stop("Group labels must not be missing.", call. = FALSE)
    levels <- unique(groups)
    references <- stats::setNames(
      lapply(levels, function(level) as.numeric(groups == level)), levels)
  }
  masses <- .massValues(.featureMetadata(data, SM.assay), rownames(metabolites))
  results <- lapply(names(references), function(reference) {
    tables <- lapply(names(matrices), function(modality) {
      expression <- matrices[[modality]]
      correlations <- vapply(seq_len(nrow(expression)), function(index) {
        values <- as.numeric(expression[index, ])
        target <- references[[reference]]
        if (length(values) < 3L || any(!is.finite(c(values, target))) ||
            stats::sd(values) == 0 || stats::sd(target) == 0) return(NA_real_)
        stats::cor(values, target, method = "pearson")
      }, numeric(1))
      data.frame(features = rownames(expression), correlation = correlations,
                 modality = modality,
                 mz = if (modality == "metabolite") masses else NA_real_)
    })
    result <- do.call(rbind, tables)
    result <- result[order(-abs(result$correlation), na.last = TRUE), ]
    result$rank <- seq_len(nrow(result))
    if (length(references) > 1L) result$ident <- reference
    if (!is.null(nfeatures)) result <- utils::head(result, nfeatures)
    result
  })
  result <- do.call(rbind, results)
  rownames(result) <- NULL
  result
}

.validateFeatureCount <- function(value, argument, allowNull = FALSE) {
  if (is.null(value) && allowNull) return(invisible(NULL))
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 1 || value > .Machine$integer.max || value != floor(value))
    stop(argument, " must be a positive integer.", call. = FALSE)
  invisible(NULL)
}

#' Find spatially variable metabolites using Moran's I
#'
#' Calls ape::Moran.I with inverse squared Euclidean distance weights and zero
#' diagonal, row-standardized by ape. Samples are analysed separately.
#' Constant features return NA.
#' @param object A SpatialExperiment or aligned MSImagingExperiment.
#' @param assay Assay name or primary/alternative experiment to analyse.
#' @param slot Assay containing expression values.
#' @param image Retained for compatibility; coordinates use spatialCoords.
#' @param nfeatures Maximum variable features selected.
#' @param max_spots Maximum pixels for the dense distance matrix; NULL uses all.
#' @param seed Random seed for sampling; caller RNG state is restored.
#' @param verbose Show progress messages.
#' @param sampleId A sample_id in colData. Required for multiple samples.
#' @return A SpatialExperiment with statistics, BH-adjusted p-values, flags and
#'   ranks in rowData of the selected experiment. Settings and sampled pixel
#'   IDs are stored in metadata(object)$moransi.
#' @export
#' @examples
#' x <- SpatialExperiment::SpatialExperiment(
#'   assays = list(counts = rbind(a = 1:6, b = c(1, 3, 5, 2, 4, 6))),
#'   rowData = S4Vectors::DataFrame(mz = c(100, 200)),
#'   spatialCoords = cbind(x = 1:6, y = 0),
#'   colData = S4Vectors::DataFrame(row.names = paste0("p", 1:6)))
#' x <- findSpatiallyVariableMetabolites(x)
#' getSpatiallyVariableMetabolites(x)
findSpatiallyVariableMetabolites <- function(
    object, assay = "main", slot = "counts", image = NULL, nfeatures = 2000,
    max_spots = 5000, seed = 1, verbose = TRUE, sampleId = NULL
) {
  if (!is.null(max_spots) &&
      (!is.numeric(max_spots) || length(max_spots) != 1L ||
       !is.finite(max_spots) || max_spots < 2 || max_spots != floor(max_spots)))
    stop("max_spots must be NULL or one finite number >= 2.", call. = FALSE)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      abs(seed) > .Machine$integer.max)
    stop("seed must be one finite number.", call. = FALSE)
  .validateFeatureCount(nfeatures, "nfeatures")
  object <- .nativeSpatialObject(object)
  expression <- .assayData(object, assay, slot)
  coordinates <- SpatialExperiment::spatialCoords(object)
  samples <- as.character(SummarizedExperiment::colData(object)$sample_id)
  if (is.null(sampleId)) {
    if (length(unique(samples)) > 1L)
      stop("Choose sampleId when analysing multiple tissue samples.", call. = FALSE)
    sampleId <- unique(samples)
  }
  if (length(sampleId) != 1L || !sampleId %in% samples)
    stop("sampleId must identify one sample in colData.", call. = FALSE)
  retained <- which(samples == sampleId)
  if (!is.null(max_spots) && length(retained) > max_spots)
    retained <- withr::with_seed(as.integer(seed), sort(sample(retained, max_spots)))
  if (length(retained) < 4L)
    stop("Moran's I requires at least four selected pixels.", call. = FALSE)
  coordinates <- coordinates[retained, , drop = FALSE]
  if (ncol(coordinates) < 2L || any(!is.finite(coordinates)))
    stop("Spatial coordinates must be finite and have at least two axes.", call. = FALSE)
  distances <- as.matrix(stats::dist(coordinates))
  if (any(distances[upper.tri(distances)] == 0))
    stop("Pixels within a sample must have distinct coordinates.", call. = FALSE)
  weights <- 1 / distances^2
  diag(weights) <- 0
  verbose_message("Computing Moran's I with inverse-square distance weights.",
                  verbose = verbose)
  statistics <- vapply(seq_len(nrow(expression)), function(index) {
    values <- as.numeric(expression[index, retained])
    if (any(!is.finite(values)) || stats::sd(values) == 0)
      return(c(observed = NA_real_, p.value = NA_real_))
    result <- ape::Moran.I(values, weights)
    c(observed = result$observed, p.value = result$p.value)
  }, c(observed = 0, p.value = 0))
  featureData <- .featureMetadata(object, assay)
  featureData$MoransI_observed <- statistics["observed", ]
  featureData$MoransI_p.value <- statistics["p.value", ]
  featureData$MoransI_p.adjust <- stats::p.adjust(statistics["p.value", ], "BH")
  valid <- which(is.finite(statistics["observed", ]) &
                   is.finite(statistics["p.value", ]))
  ordered <- valid[order(statistics["p.value", valid],
                         -abs(statistics["observed", valid]))]
  ranks <- rep(NA_integer_, nrow(expression))
  ranks[ordered] <- seq_along(ordered)
  featureData$moransi.spatially.variable.rank <- ranks
  featureData$moransi.spatially.variable <- !is.na(ranks) & ranks <= nfeatures
  object <- .setFeatureMetadata(object, featureData, assay)
  S4Vectors::metadata(object)$moransi <- list(
    backend = "ape::Moran.I", weights = "row-standardized inverse squared Euclidean distance",
    sample_id = sampleId, assay = assay, layer = slot,
    pixels = colnames(object)[retained], seed = seed)
  object
}

#' Get the top spatially variable metabolites
#' @param object A SpatialExperiment returned by findSpatiallyVariableMetabolites.
#' @param assay Primary or alternative experiment containing results.
#' @param n Maximum number of feature names.
#' @return Feature names in Moran's I rank order, excluding untested features.
#' @export
getSpatiallyVariableMetabolites <- function(object, assay = "main", n = 10) {
  .validateFeatureCount(n, "n")
  metadata <- .featureMetadata(object, assay)
  ranks <- metadata$moransi.spatially.variable.rank
  if (is.null(ranks))
    stop("Run findSpatiallyVariableMetabolites() first.", call. = FALSE)
  tested <- which(!is.na(ranks))
  utils::head(rownames(metadata)[tested[order(ranks[tested])]], n)
}

list_to_pprcomp <- function(lst) {
  structure(lst[c("sdev", "rotation", "center", "scale", "x")], class = "prcomp")
}

RowVar <- function(x) {
  if (ncol(x) < 2L) return(rep(NA_real_, nrow(x)))
  means <- rowMeans(x)
  pmax(0, (rowSums(x * x) - ncol(x) * means * means) / (ncol(x) - 1L))
}
