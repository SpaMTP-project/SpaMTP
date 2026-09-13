#' Mult-Omic data integration
#'
#' SpatialExperiment input uses paired alternative experiments and combines
#' equal-weight PCA embeddings generated with scater. This is not Seurat WNN;
#' the historical WNN workflow remains in the published-workflow branch.
#'
#' @param multiomic.data SpaMTP dataset contain Spatial Transcriptomics and Metabolomic datasets in two different assays
#' @param reduction.list Reduction names for the primary MSI and alternative
#'   modalities (default = list("spm.pca", "spt.pca")).
#' @param dims.list List containing the numeric range of principle component dimension to include for each modality (default = list(1:30,1:30)).
#' @param return.intermediate Retain per-modality PCA results in reducedDims().
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#' @param modalities Primary experiment (`"main"`) followed by alternative
#'   experiments to integrate. Derived pathway and merged assays are not
#'   included automatically.
#' @param ... Additional arguments passed to scater::runPCA().
#'
#' @return The input container with an integrated representation. For
#'   SpatialExperiment this is `reducedDim(x, "integrated")`.
#' @export
#'
#' @examples
#' utils::str(formals(multiOmicIntegration))
#' # SpaMTP.obj <- multiOmicIntegration(SpaMTP.obj, reduction.list =  list("spt.pca", "spm.pca"), dims.list = list(1:30, 1:30))
methods::setGeneric(
  "multiOmicIntegration",
  function(
      multiomic.data,
      reduction.list = list("spm.pca", "spt.pca"),
      dims.list = list(1:30, 1:30),
      return.intermediate = FALSE,
      verbose = FALSE,
      modalities = c("main", "transcriptome"),
      ...
  ) {
    methods::standardGeneric("multiOmicIntegration")
  }
)

.integrationAssay <- function(experiment) {
  available <- SummarizedExperiment::assayNames(experiment)
  selected <- intersect(c("logcounts", "normcounts", "counts"), available)
  if (!length(selected)) {
    stop(
      "Each modality needs one of these assays: logcounts, normcounts, counts.",
      call. = FALSE
    )
  }
  selected[[1L]]
}

.modalityPca <- function(experiment, reductionName, dimensions, ...) {
  if (!length(dimensions) || anyNA(dimensions) || any(dimensions < 1L) ||
      anyDuplicated(dimensions)) {
    stop("Each dims.list entry must contain unique positive dimensions.", call. = FALSE)
  }
  if (!methods::is(experiment, "SingleCellExperiment")) {
    experiment <- methods::as(experiment, "SingleCellExperiment")
  }
  available <- SingleCellExperiment::reducedDimNames(experiment)
  if (!reductionName %in% available) {
    maximum <- min(nrow(experiment), ncol(experiment)) - 1L
    if (maximum < 1L) {
      stop("Each modality needs at least two features and pixels.", call. = FALSE)
    }
    ncomponents <- min(max(as.integer(dimensions)), maximum)
    experiment <- scater::runPCA(
      experiment,
      exprs_values = .integrationAssay(experiment),
      ncomponents = ncomponents,
      name = reductionName,
      ...
    )
  }
  embedding <- SingleCellExperiment::reducedDim(experiment, reductionName)
  retained <- dimensions[dimensions <= ncol(embedding)]
  if (!length(retained)) {
    stop("Requested dimensions are absent from reduction `", reductionName, "`.")
  }
  list(experiment = experiment, embedding = embedding[, retained, drop = FALSE])
}

#' @rdname multiOmicIntegration
#' @export
methods::setMethod(
  "multiOmicIntegration",
  "SpatialExperiment",
  function(
      multiomic.data,
      reduction.list = list("spm.pca", "spt.pca"),
      dims.list = list(1:30, 1:30),
      return.intermediate = FALSE,
      verbose = FALSE,
      modalities = c("main", "transcriptome"),
      ...
  ) {
    if (length(modalities) < 2L || !identical(modalities[[1L]], "main") ||
        anyNA(modalities) || anyDuplicated(modalities)) {
      stop("modalities must start with main and contain distinct altExp names.",
           call. = FALSE)
    }
    alternatives <- modalities[-1L]
    if (!all(alternatives %in% SingleCellExperiment::altExpNames(multiomic.data))) {
      stop(
        "Requested alternative modality was not found. Use addTranscriptome() ",
        "or supply modalities matching altExpNames(x).",
        call. = FALSE
      )
    }
    modalityNames <- c("main", alternatives)
    experiments <- c(
      list(main = multiomic.data),
      stats::setNames(
        lapply(alternatives, function(name) {
          SingleCellExperiment::altExp(multiomic.data, name)
        }),
        alternatives
      )
    )
    if (length(reduction.list) != length(experiments)) {
      stop("Provide one reduction.list entry per modality.", call. = FALSE)
    }
    if (length(dims.list) != length(experiments)) {
      stop("Provide one dims.list entry per modality.", call. = FALSE)
    }

    results <- lapply(seq_along(experiments), function(index) {
      .modalityPca(
        experiments[[index]],
        reductionName = as.character(reduction.list[[index]])[[1L]],
        dimensions = as.integer(dims.list[[index]]),
        ...
      )
    })
    embeddings <- lapply(seq_along(results), function(index) {
      embedding <- scale(results[[index]]$embedding)
      embedding[!is.finite(embedding)] <- 0
      embedding <- embedding / sqrt(ncol(embedding))
      colnames(embedding) <- paste(
        modalityNames[[index]],
        colnames(embedding) %||% seq_len(ncol(embedding)),
        sep = "_"
      )
      embedding
    })
    integrated <- do.call(cbind, embeddings)
    SingleCellExperiment::reducedDim(multiomic.data, "integrated") <- integrated

    if (isTRUE(return.intermediate)) {
      mainReductions <- SingleCellExperiment::reducedDims(results[[1L]]$experiment)
      SingleCellExperiment::reducedDims(multiomic.data) <- mainReductions
      SingleCellExperiment::reducedDim(multiomic.data, "integrated") <- integrated
      for (index in seq_along(alternatives)) {
        SingleCellExperiment::altExp(multiomic.data, alternatives[[index]]) <-
          results[[index + 1L]]$experiment
      }
    }
    S4Vectors::metadata(multiomic.data)$spamtp_integration <- list(
      method = "equal-weight concatenated PCA",
      modalities = modalityNames,
      reductions = unlist(reduction.list),
      dimensions = dims.list
    )
    verbose_message(
      paste0(
        "Stored a ", ncol(integrated),
        "-dimensional joint embedding in reducedDim(x, `integrated`)."
      ),
      verbose = verbose
    )
    multiomic.data
  }
)

#' @rdname multiOmicIntegration
#' @export
methods::setMethod(
  "multiOmicIntegration",
  "ANY",
  function(
      multiomic.data,
      reduction.list = list("spm.pca", "spt.pca"),
      dims.list = list(1:30, 1:30),
      return.intermediate = FALSE,
      verbose = FALSE,
      modalities = c("main", "transcriptome"),
      ...
  ) {
    .requireExperiment(multiomic.data, "SpatialExperiment")
  }
)




#' Create a singular multiomics assay by merging data from multiple assays.
#'
#' Combines selected modalities as a scaled alternative experiment.
#' Useful for integrating multiple modalities (e.g. transcriptomics, proteomics, metabolomics) that have already been scaled.
#'
#' @param SpaMTP A Bioconductor experiment that contains atleast two assays to be merged.
#' @param assays.to.merge At least two distinct modality names: main and/or altExp names.
#' @param new.assay A character string specifying the name of the new assay to be created (default = "merged").
#' @param return.original TRUE adds an altExp; FALSE returns only the merged feature space.
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return A Bioconductor experiment containing a new assay with the merged scaled data values.
#' @export
#'
#' @details
#' Uses an existing scaled assay when available; otherwise centres and scales
#' logcounts, normcounts or counts across pixels. Constant features become zero.
#' The merged values are stored as scaled, never relabelled as counts.
#'
#' @examples
#' utils::str(formals(createMergedModalityAssay))
#' # merged_obj <- createMergedModalityAssay(SpaMTP = spamtp_obj, assays.to.merge = c("SPM", "SPT"),new.assay = "merged")
createMergedModalityAssay <- function(SpaMTP, assays.to.merge, new.assay = "merged", return.original = TRUE, verbose = FALSE){
  .requireExperiment(SpaMTP, "SpatialExperiment")
  modalityNames <- SingleCellExperiment::altExpNames(SpaMTP)
  primaryNames <- unique(c(
    "main", "primary", SummarizedExperiment::assayNames(SpaMTP)
  ))
  getModality <- function(name) {
    if (name %in% modalityNames) {
      return(SingleCellExperiment::altExp(SpaMTP, name))
    }
    if (name %in% primaryNames) {
      return(SpaMTP)
    }
    stop(
      "Modality `", name, "` was not found. Use `main` or an altExp name.",
      call. = FALSE
    )
  }
  if (length(assays.to.merge) < 2L || anyNA(assays.to.merge) ||
      anyDuplicated(assays.to.merge)) {
    stop("At least two distinct modalities must be supplied.", call. = FALSE)
  }
  modalities <- lapply(assays.to.merge, getModality)
  matrices <- Map(
    function(modality, name) {
      available <- SummarizedExperiment::assayNames(modality)
      preferred <- c("scale.data", "scaled", "logcounts", "normcounts", "counts")
      selected <- intersect(preferred, available)
      if (!length(selected)) stop("No supported expression assay in modality ", name, call. = FALSE)
      selected <- selected[[1L]]
      matrix <- SummarizedExperiment::assay(modality, selected)
      if (!selected %in% c("scale.data", "scaled")) {
        matrix <- t(scale(t(as.matrix(matrix))))
        matrix[!is.finite(matrix)] <- 0
      }
      rownames(matrix) <- make.unique(paste(name, rownames(matrix), sep = "::"))
      matrix
    },
    modalities,
    assays.to.merge
  )
  merged <- do.call(rbind, matrices)
  featureData <- S4Vectors::DataFrame(
    modality = rep(assays.to.merge, vapply(matrices, nrow, integer(1))),
    original_feature = unlist(lapply(matrices, function(x) sub("^[^:]+::", "", rownames(x))))
  )
  rownames(featureData) <- rownames(merged)
  mergedExperiment <- SingleCellExperiment::SingleCellExperiment(
    assays = list(scaled = merged),
    rowData = featureData
  )
  if (isTRUE(return.original)) {
    SingleCellExperiment::altExp(SpaMTP, new.assay) <- mergedExperiment
    return(SpaMTP)
  }
  return(SpatialExperiment::SpatialExperiment(
    assays = list(scaled = merged),
    rowData = featureData,
    colData = SummarizedExperiment::colData(SpaMTP),
    spatialCoords = SpatialExperiment::spatialCoords(SpaMTP),
    imgData = SpatialExperiment::imgData(SpaMTP),
    metadata = S4Vectors::metadata(SpaMTP)
  ))
}
