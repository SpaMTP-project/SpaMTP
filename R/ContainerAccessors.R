.requireExperiment <- function(object, class = "SummarizedExperiment") {
  if (!methods::is(object, class)) {
    stop(
      "Use a ", class, ". Convert Seurat input explicitly with ",
      "seuratToSpatialExperiment() or seuratToSingleCellExperiment() first.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.isSummarizedExperiment <- function(object) {
  inherits(object, "SummarizedExperiment")
}

.experimentForAssay <- function(object, assay = NULL) {
  .requireExperiment(object)
  alternatives <- if (inherits(object, "SingleCellExperiment")) {
    SingleCellExperiment::altExpNames(object)
  } else {
    character()
  }
  if (!is.null(assay) && assay %in% alternatives) {
    return(SingleCellExperiment::altExp(object, assay))
  }
  primary <- c("main", "primary", "Spatial", "SPM",
               SummarizedExperiment::assayNames(object))
  if (inherits(object, "SingleCellExperiment")) {
    primary <- c(primary, SingleCellExperiment::mainExpName(object))
  }
  if (!is.null(assay) && !assay %in% primary) {
    stop(
      "Unknown assay or modality `", assay, "`. Use `main`, an assay name, ",
      "or an altExp name: ", paste(alternatives, collapse = ", "), ".",
      call. = FALSE
    )
  }
  object
}

.replaceExperiment <- function(object, experiment, assay = NULL) {
  if (!is.null(assay) && inherits(object, "SingleCellExperiment") &&
      assay %in% SingleCellExperiment::altExpNames(object)) {
    SingleCellExperiment::altExp(object, assay) <- experiment
    return(object)
  }
  experiment
}

.alignMetadata <- function(value, target, what) {
  value <- S4Vectors::DataFrame(value)
  identifiers <- rownames(value)
  if (is.null(identifiers) ||
      identical(identifiers, as.character(seq_len(nrow(value))))) {
    if (nrow(value) != length(target)) {
      stop(what, " must contain one row for every identifier.", call. = FALSE)
    }
    rownames(value) <- target
  } else if (anyDuplicated(identifiers) ||
             !setequal(identifiers, target)) {
    stop(what, " row names must match the object identifiers.", call. = FALSE)
  }
  value[target, , drop = FALSE]
}

.nativeSpatialObject <- function(object) {
  if (methods::is(object, "MSImagingExperiment")) {
    return(asSpatialExperiment(object))
  }
  if (!methods::is(object, "SpatialExperiment")) {
    stop(
      "Use a SpatialExperiment or aligned MSImagingExperiment; convert ",
      "Seurat input explicitly with seuratToSpatialExperiment().",
      call. = FALSE
    )
  }
  object
}

.cellMetadata <- function(object) {
  .requireExperiment(object)
  as.data.frame(SummarizedExperiment::colData(object), optional = TRUE)
}

.setCellMetadata <- function(object, value) {
  .requireExperiment(object)
  value <- .alignMetadata(value, colnames(object), "Cell metadata")
  existing <- SummarizedExperiment::colData(object)
  for (column in colnames(value)) {
    existing[[column]] <- value[[column]]
  }
  SummarizedExperiment::colData(object) <- existing
  object
}

.featureMetadata <- function(object, assay = NULL) {
  experiment <- .experimentForAssay(object, assay)
  as.data.frame(SummarizedExperiment::rowData(experiment), optional = TRUE)
}

.setFeatureMetadata <- function(object, value, assay = NULL) {
  experiment <- .experimentForAssay(object, assay)
  SummarizedExperiment::rowData(experiment) <-
    .alignMetadata(value, rownames(experiment), "Feature metadata")
  .replaceExperiment(object, experiment, assay)
}

.assayData <- function(object, assay = NULL, layer = "counts") {
  experiment <- .experimentForAssay(object, assay)
  available <- SummarizedExperiment::assayNames(experiment)
  if (length(layer) != 1L || is.na(layer) || !layer %in% available) {
    stop(
      "Assay `", paste(layer, collapse = ", "), "` was not found. Available assays: ",
      paste(available, collapse = ", "), ".", call. = FALSE
    )
  }
  SummarizedExperiment::assay(experiment, layer)
}

.setAssayData <- function(object, value, assay = NULL, layer = "counts") {
  experiment <- .experimentForAssay(object, assay)
  if (length(layer) != 1L || is.na(layer) || !nzchar(layer)) {
    stop("layer must be one non-empty assay name.", call. = FALSE)
  }
  SummarizedExperiment::assay(experiment, layer) <- value
  .replaceExperiment(object, experiment, assay)
}

.assayNames <- function(object) {
  .requireExperiment(object)
  SummarizedExperiment::assayNames(object)
}

.storedDataNames <- function(object) {
  .requireExperiment(object)
  names(S4Vectors::metadata(object))
}

.storedData <- function(object, name) {
  .requireExperiment(object)
  S4Vectors::metadata(object)[[name]]
}

.setStoredData <- function(object, name, value) {
  .requireExperiment(object)
  metadata <- S4Vectors::metadata(object)
  metadata[[name]] <- value
  S4Vectors::metadata(object) <- metadata
  object
}
