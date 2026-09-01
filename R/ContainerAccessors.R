.isSeuratObject <- function(object) {
  inherits(object, "Seurat")
}

.isSummarizedExperiment <- function(object) {
  inherits(object, "SummarizedExperiment")
}

.cellMetadata <- function(object) {
  if (.isSeuratObject(object)) {
    return(object[[]])
  }
  if (.isSummarizedExperiment(object)) {
    return(as.data.frame(SummarizedExperiment::colData(object)))
  }
  stop("Unsupported container class: ", paste(class(object), collapse = "/"), call. = FALSE)
}

.setCellMetadata <- function(object, value) {
  value <- as.data.frame(value)
  target <- colnames(object)
  valueRows <- rownames(value)
  defaultRows <- identical(valueRows, as.character(seq_len(nrow(value))))
  if (is.null(valueRows) || defaultRows || !all(target %in% valueRows)) {
    if (nrow(value) != length(target)) {
      stop(
        "Cell metadata must contain one row for every cell in the object.",
        call. = FALSE
      )
    }
    rownames(value) <- target
  }
  value <- value[target, , drop = FALSE]
  if (.isSeuratObject(object)) {
    return(SeuratObject::AddMetaData(object, metadata = value))
  }
  if (.isSummarizedExperiment(object)) {
    SummarizedExperiment::colData(object) <- S4Vectors::DataFrame(value)
    return(object)
  }
  stop("Unsupported container class: ", paste(class(object), collapse = "/"), call. = FALSE)
}

.featureMetadata <- function(object, assay = NULL) {
  if (.isSeuratObject(object)) {
    assay <- assay %||% SeuratObject::DefaultAssay(object)
    return(object[[assay]][[]])
  }
  if (.isSummarizedExperiment(object)) {
    return(as.data.frame(SummarizedExperiment::rowData(object)))
  }
  stop("Unsupported container class: ", paste(class(object), collapse = "/"), call. = FALSE)
}

.setFeatureMetadata <- function(object, value, assay = NULL) {
  value <- as.data.frame(value)
  if (.isSeuratObject(object)) {
    assay <- assay %||% SeuratObject::DefaultAssay(object)
    assayObject <- object[[assay]]
    target <- rownames(assayObject)
    valueRows <- rownames(value)
    defaultRows <- identical(valueRows, as.character(seq_len(nrow(value))))
    if (is.null(valueRows) || defaultRows || !all(target %in% valueRows)) {
      if (nrow(value) != length(target)) {
        stop(
          "Feature metadata must contain one row for every feature in the assay.",
          call. = FALSE
        )
      }
      rownames(value) <- target
    }
    value <- value[target, , drop = FALSE]
    for (column in setdiff(colnames(assayObject[[]]), colnames(value))) {
      assayObject[[column]] <- NULL
    }
    if (ncol(value)) {
      assayObject[[]] <- value
    }
    object[[assay]] <- assayObject
    return(object)
  }
  if (.isSummarizedExperiment(object)) {
    target <- rownames(object)
    valueRows <- rownames(value)
    defaultRows <- identical(valueRows, as.character(seq_len(nrow(value))))
    if (is.null(valueRows) || defaultRows || !all(target %in% valueRows)) {
      if (nrow(value) != length(target)) {
        stop(
          "Feature metadata must contain one row for every feature in the object.",
          call. = FALSE
        )
      }
      rownames(value) <- target
    }
    value <- value[target, , drop = FALSE]
    SummarizedExperiment::rowData(object) <- S4Vectors::DataFrame(value)
    return(object)
  }
  stop("Unsupported container class: ", paste(class(object), collapse = "/"), call. = FALSE)
}

.assayData <- function(object, assay = NULL, layer = "counts") {
  if (.isSeuratObject(object)) {
    assay <- assay %||% SeuratObject::DefaultAssay(object)
    return(SeuratObject::LayerData(object, assay = assay, layer = layer))
  }
  if (.isSummarizedExperiment(object)) {
    available <- SummarizedExperiment::assayNames(object)
    selected <- if (layer %in% available) layer else assay
    if (is.null(selected) || !selected %in% available) {
      stop(
        "Assay `", layer, "` was not found. Available assays: ",
        paste(available, collapse = ", "), ".",
        call. = FALSE
      )
    }
    return(SummarizedExperiment::assay(object, selected))
  }
  stop("Unsupported container class: ", paste(class(object), collapse = "/"), call. = FALSE)
}

.setAssayData <- function(object, value, assay = NULL, layer = "counts") {
  if (.isSeuratObject(object)) {
    assay <- assay %||% SeuratObject::DefaultAssay(object)
    SeuratObject::LayerData(object, assay = assay, layer = layer) <- value
    return(object)
  }
  if (.isSummarizedExperiment(object)) {
    available <- SummarizedExperiment::assayNames(object)
    selected <- if (layer %in% available) layer else assay
    if (is.null(selected) || !selected %in% available) {
      stop(
        "Assay `", layer, "` was not found. Available assays: ",
        paste(available, collapse = ", "), ".",
        call. = FALSE
      )
    }
    SummarizedExperiment::assay(object, selected) <- value
    return(object)
  }
  stop("Unsupported container class: ", paste(class(object), collapse = "/"), call. = FALSE)
}

.assayNames <- function(object) {
  if (.isSeuratObject(object)) {
    return(SeuratObject::Assays(object))
  }
  if (.isSummarizedExperiment(object)) {
    return(SummarizedExperiment::assayNames(object))
  }
  stop("Unsupported container class: ", paste(class(object), collapse = "/"), call. = FALSE)
}

.storedDataNames <- function(object) {
  if (!.isSeuratObject(object)) {
    return(names(S4Vectors::metadata(object)))
  }
  unique(c(
    names(SeuratObject::Misc(object)),
    SeuratObject::Tool(object)
  ))
}

.storedData <- function(object, name) {
  if (!.isSeuratObject(object)) {
    return(S4Vectors::metadata(object)[[name]])
  }
  misc <- SeuratObject::Misc(object)
  if (name %in% names(misc)) {
    value <- misc[[name]]
    if (is.list(value) && identical(names(value), ".spamtp_data_frame")) {
      return(value[[1L]])
    }
    return(value)
  }
  tryCatch(
    SeuratObject::Tool(object, slot = name),
    error = function(e) NULL
  )
}

.setStoredData <- function(object, name, value) {
  if (!.isSeuratObject(object)) {
    metadata <- S4Vectors::metadata(object)
    metadata[[name]] <- value
    S4Vectors::metadata(object) <- metadata
    return(object)
  }
  if (is.data.frame(value)) {
    value <- list(.spamtp_data_frame = value)
  }
  SeuratObject::Misc(object, slot = name) <- value
  object
}

.copyStoredData <- function(from, to) {
  for (name in .storedDataNames(from)) {
    to <- .setStoredData(to, name, .storedData(from, name))
  }
  to
}

.centroidsWithCoordinates <- function(centroids, coordinates) {
  cells <- as.character(SeuratObject::Cells(centroids))
  coordinates <- as.data.frame(coordinates)
  if (!"cell" %in% names(coordinates)) {
    coordinates$cell <- cells
  }
  coordinates <- coordinates[match(cells, coordinates$cell), c("x", "y", "cell"), drop = FALSE]
  if (anyNA(coordinates[c("x", "y", "cell")])) {
    stop("Coordinates must contain one x/y pair for every centroid.", call. = FALSE)
  }
  SeuratObject::CreateCentroids(
    coordinates,
    nsides = length(centroids),
    radius = SeuratObject::Radius(centroids),
    theta = SeuratObject::Theta(centroids)
  )
}

.fovWithBoundary <- function(fov, boundary, value) {
  boundaryNames <- SeuratObject::Boundaries(fov)
  boundaries <- stats::setNames(
    lapply(boundaryNames, function(name) fov[[name]]),
    boundaryNames
  )
  boundaries[[boundary]] <- value
  updated <- SeuratObject::CreateFOV(
    coords = boundaries,
    molecules = NULL,
    assay = SeuratObject::DefaultAssay(fov),
    key = SeuratObject::Key(fov)
  )
  for (name in SeuratObject::Molecules(fov)) {
    updated[[name]] <- fov[[name]]
  }
  updated
}
