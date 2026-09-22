.seuratLayer <- function(x, assay, layer) {
  if (length(assay) != 1L || is.na(assay) ||
      !assay %in% SeuratObject::Assays(x)) {
    stop("Select an existing Seurat assay.", call. = FALSE)
  }
  available <- SeuratObject::Layers(x[[assay]], search = NA)
  if (length(layer) != 1L || is.na(layer) || !layer %in% available) {
    stop(
      "Select one exact Seurat layer. Available layers: ",
      paste(available, collapse = ", "),
      ". Join split layers in Seurat first, or select one explicitly.",
      call. = FALSE
    )
  }
  SeuratObject::LayerData(x, assay = assay, layer = layer)
}

.hasSeuratMetadata <- function(value) {
  if (inherits(value, "Seurat") ||
      any(attr(class(value), "package") %in% c("Seurat", "SeuratObject"))) {
    return(TRUE)
  }
  is.list(value) && any(vapply(value, .hasSeuratMetadata, logical(1)))
}

.seuratStoredMetadata <- function(x) {
  result <- SeuratObject::Misc(x)
  for (name in setdiff(SeuratObject::Tool(x), names(result))) {
    result[[name]] <- SeuratObject::Tool(x, slot = name)
  }
  result <- lapply(result, function(value) {
    if (is.list(value) && identical(names(value), ".spamtp_data_frame")) {
      value[[1L]]
    } else value
  })
  unsupported <- vapply(result, .hasSeuratMetadata, logical(1))
  if (any(unsupported)) {
    warning("Skipping Seurat-specific stored metadata: ",
            paste(names(result)[unsupported], collapse = ", "), call. = FALSE)
  }
  result[!unsupported]
}

.seuratAssayToExperiment <- function(x, assay, layer = "counts") {
  matrix <- .seuratLayer(x, assay, layer)
  featureData <- x[[assay]][[]]
  if (is.null(rownames(matrix)) || is.null(colnames(matrix)) ||
      anyDuplicated(rownames(matrix)) || anyDuplicated(colnames(matrix))) {
    stop("Converted layers need unique feature and pixel names.", call. = FALSE)
  }
  SingleCellExperiment::SingleCellExperiment(
    assays = stats::setNames(list(matrix), layer),
    rowData = S4Vectors::DataFrame(
      featureData[rownames(matrix), , drop = FALSE], check.names = FALSE),
    colData = S4Vectors::DataFrame(
      x[[]][colnames(matrix), , drop = FALSE], check.names = FALSE)
  )
}

#' Convert Seurat expression data to SingleCellExperiment
#'
#' Conversion is the only part of SpaMTP that uses SeuratObject. Coordinates
#' are not required. The selected layer keeps its original name (for example,
#' `counts` or `data`); select that assay explicitly in downstream functions.
#' Other Seurat assays with the same pixels become `altExp()` entries.
#' Only one exact layer is converted, not all split layers. Join split layers
#' in Seurat before conversion when the complete dataset is needed.
#' Cell/feature metadata, active identities and stored analysis metadata are
#' retained. Reductions, graphs, optical images and Seurat commands are not
#' converted. Stored metadata should contain portable R/Bioconductor values;
#' entries containing Seurat objects/classes in lists are skipped with a warning.
#'
#' @param x A Seurat object.
#' @param assay Primary Seurat assay; defaults to its active assay.
#' @param layer Exact layer to convert. No partial matching is performed.
#' @param includeAltExps Convert other assays with the same pixels.
#'
#' @return A SingleCellExperiment with the selected expression layer.
#' @export
#' @examples
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'     counts <- matrix(seq_len(12), nrow = 3,
#'         dimnames = list(paste0("gene", 1:3), paste0("cell", 1:4)))
#'     seu <- SeuratObject::CreateSeuratObject(counts)
#'     sce <- seuratToSingleCellExperiment(seu)
#'     SummarizedExperiment::assayNames(sce)
#' }
seuratToSingleCellExperiment <- function(
    x, assay = NULL, layer = "counts", includeAltExps = TRUE
) {
  .requireOptionalPackage("SeuratObject", "convert a Seurat object")
  if (!inherits(x, "Seurat")) {
    stop("`x` must be a Seurat object.", call. = FALSE)
  }
  assay <- assay %||% SeuratObject::DefaultAssay(x)
  object <- .seuratAssayToExperiment(x, assay, layer)
  SingleCellExperiment::mainExpName(object) <- assay
  labels <- SeuratObject::Idents(x)[colnames(object)]
  SingleCellExperiment::colLabels(object) <- labels
  S4Vectors::metadata(object) <- .seuratStoredMetadata(x)
  if (isTRUE(includeAltExps)) {
    for (alternative in setdiff(SeuratObject::Assays(x), assay)) {
      experiment <- tryCatch(
        .seuratAssayToExperiment(x, alternative, layer),
        error = function(error) {
          warning("Skipping Seurat assay ", alternative, ": ",
                  conditionMessage(error), call. = FALSE)
          NULL
        }
      )
      if (is.null(experiment)) next
      order <- match(colnames(object), colnames(experiment))
      if (!anyNA(order) && ncol(experiment) == ncol(object)) {
        SingleCellExperiment::altExp(object, alternative) <-
          experiment[, order, drop = FALSE]
      } else {
        warning("Skipping unpaired Seurat assay: ", alternative, call. = FALSE)
      }
    }
  }
  S4Vectors::metadata(object)$seurat_interoperability <- list(
    primary_assay = assay, layer = layer)
  object
}

.seuratCoordinates <- function(x, image = NULL, cells = colnames(x),
                               coordinates = NULL) {
  if (!is.null(coordinates) && !is.null(image)) {
    stop("Supply coordinates or image, not both.", call. = FALSE)
  }
  source <- "coordinates"
  columns <- c("x", "y")
  if (is.null(coordinates) && is.null(image)) {
    images <- SeuratObject::Images(x)
    if (length(images) > 1L) {
      stop("Multiple Seurat images found. Supply image or named coordinates.",
           call. = FALSE)
    }
    if (length(images) == 1L) {
      # FOV coordinates may have been aligned after the metadata was created.
      # Never silently replace the spatial accessor's values with that cache.
      image <- images[[1L]]
    } else {
      pixelData <- x[[]]
      columns <- if (all(c("x", "y") %in% names(pixelData))) c("x", "y") else
        if (all(c("x_coord", "y_coord") %in% names(pixelData))) c("x_coord", "y_coord") else NULL
      if (!is.null(columns)) {
        coordinates <- pixelData[, columns, drop = FALSE]
        colnames(coordinates) <- c("x", "y")
        source <- "cell_metadata"
      }
    }
  }
  if (!is.null(image)) {
    if (length(image) != 1L || is.na(image) ||
        !image %in% SeuratObject::Images(x)) {
      stop("Select an existing Seurat image.", call. = FALSE)
    }
    coordinates <- tryCatch(
      SeuratObject::GetTissueCoordinates(x, image = image, full = FALSE),
      error = function(error) {
        stop("Cannot read image coordinates with GetTissueCoordinates(): ",
             conditionMessage(error), ". Update the Seurat object or supply ",
             "named x/y coordinates explicitly.", call. = FALSE)
      })
    source <- "image"
  }
  if (is.null(coordinates)) {
    stop(
      "No spatial coordinates found. Supply named x/y coordinates, or use ",
      "seuratToSingleCellExperiment() for non-spatial data.", call. = FALSE)
  }
  coordinates <- as.data.frame(coordinates)
  if (!all(c("x", "y") %in% names(coordinates))) {
    stop("Coordinates must contain x and y columns; map image axes explicitly.",
         call. = FALSE)
  }
  identifiers <- if ("cell" %in% names(coordinates)) {
    as.character(coordinates$cell)
  } else rownames(coordinates)
  order <- match(cells, identifiers)
  if (is.null(identifiers) || anyNA(identifiers) ||
      anyDuplicated(identifiers) || anyNA(order)) {
    stop("Coordinates must uniquely identify every selected pixel.", call. = FALSE)
  }
  result <- as.matrix(coordinates[order, c("x", "y"), drop = FALSE])
  if (!is.numeric(result) || any(!is.finite(result))) {
    stop("Coordinates must be finite numeric x/y values.", call. = FALSE)
  }
  rownames(result) <- cells
  list(values = result, source = source, image = image, columns = columns)
}

#' Convert a Seurat object to SpatialExperiment
#'
#' Convert once at the workflow boundary, then use Bioconductor containers
#' throughout SpaMTP. Expression and metadata follow
#' [seuratToSingleCellExperiment()]. Coordinates are matched by pixel names,
#' never by position. Supply `coordinates` for multi-image data or image
#' formats without x/y columns; coordinate units and orientation must already
#' be consistent. Optical rasters, graphs and reductions are not transferred.
#' Attach optical images separately with [addSpatialImage()].
#' When a single image/FOV is present, its public `GetTissueCoordinates()`
#' accessor is the default coordinate source, even if historical x/y columns
#' remain in cell metadata. Metadata coordinates are used automatically only
#' when there are no images. Multiple images always require an explicit
#' selection. Centroids are read as pixel centres, not expanded polygon
#' vertices. The selected source is recorded in
#' `metadata(x)$seurat_interoperability`; retained legacy coordinate columns
#' are archival and must not replace `spatialCoords(x)` in native analysis.
#'
#' @param x A Seurat object.
#' @param assay Primary Seurat assay; defaults to its active assay.
#' @param layer Exact layer to convert, usually `"counts"`.
#' @param image Optional Seurat image/FOV used for coordinates. A single image
#'   is selected automatically; its coordinates take precedence over metadata.
#' @param includeAltExps Convert other assays with matching pixels.
#' @param coordinates Optional matrix/data.frame with x/y columns and pixel
#'   names in row names or a `cell` column. Cannot be combined with image.
#'
#' @return A SpatialExperiment.
#' @export
#' @examples
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'     counts <- matrix(seq_len(12), nrow = 3,
#'         dimnames = list(paste0("mz-", 101:103), paste0("pixel", 1:4)))
#'     seu <- SeuratObject::CreateSeuratObject(counts)
#'     xy <- cbind(x = 1:4, y = c(0, 1, 0, 1))
#'     rownames(xy) <- colnames(seu)
#'     spe <- seuratToSpatialExperiment(seu, coordinates = xy)
#'     spe <- normalizeSMData(spe, verbose = FALSE)
#' }
seuratToSpatialExperiment <- function(
    x, assay = NULL, layer = "counts", image = NULL,
    includeAltExps = TRUE, coordinates = NULL
) {
  experiment <- seuratToSingleCellExperiment(x, assay, layer, includeAltExps)
  coordinateInput <- .seuratCoordinates(
    x, image, cells = colnames(experiment), coordinates = coordinates)
  pixelData <- SummarizedExperiment::colData(experiment)
  samples <- if ("sample_id" %in% colnames(pixelData)) {
    as.character(pixelData$sample_id)
  } else if ("orig.ident" %in% colnames(pixelData)) {
    as.character(pixelData$orig.ident)
  } else rep("sample01", ncol(experiment))
  object <- SpatialExperiment::SpatialExperiment(
    assays = SummarizedExperiment::assays(experiment),
    rowData = SummarizedExperiment::rowData(experiment),
    colData = pixelData, spatialCoords = coordinateInput$values, sample_id = samples,
    metadata = S4Vectors::metadata(experiment))
  SingleCellExperiment::mainExpName(object) <-
    SingleCellExperiment::mainExpName(experiment)
  for (name in SingleCellExperiment::altExpNames(experiment)) {
    SingleCellExperiment::altExp(object, name) <-
      SingleCellExperiment::altExp(experiment, name)
  }
  provenance <- S4Vectors::metadata(object)$seurat_interoperability
  provenance$image <- coordinateInput$image
  provenance$coordinate_source <- coordinateInput$source
  provenance$coordinate_columns <- coordinateInput$columns
  S4Vectors::metadata(object)$seurat_interoperability <- provenance
  object
}

#' Convert a SpatialExperiment to an optional Seurat object
#'
#' This conversion is provided only for users who explicitly need a Seurat
#' workflow. Only SeuratObject is needed for conversion; Seurat analysis is not used.
#' The requested layer is placed in the Seurat counts layer; choose counts
#' when raw intensities are needed. Other layers, reductions and optical
#' images are not converted; spatial coordinates are exported as centroids.
#' Alternative experiments without the requested layer are skipped with a
#' warning, rather than treating a derived score assay as raw counts.
#'
#' @param x A `SpatialExperiment`.
#' @param assayName Name for the primary Seurat assay.
#' @param layer SpatialExperiment assay to convert.
#' @param imageName Name for the centroid FOV.
#' @param includeAltExps Convert alternative experiments to Seurat assays.
#'
#' @return A Seurat object.
#' @export
#' @examples
#' utils::str(formals(spatialExperimentToSeurat))
spatialExperimentToSeurat <- function(
    x,
    assayName = "Spatial",
    layer = "counts",
    imageName = "fov",
    includeAltExps = TRUE
) {
  .requireOptionalPackage("SeuratObject", "convert to a Seurat object")
  if (!methods::is(x, "SpatialExperiment")) {
    stop("`x` must be a SpatialExperiment.", call. = FALSE)
  }
  matrix <- .assayData(x, layer = layer)
  if (!inherits(matrix, "sparseMatrix")) {
    matrix <- Matrix::Matrix(matrix, sparse = TRUE)
  }
  object <- SeuratObject::CreateSeuratObject(counts = matrix, assay = assayName)
  object <- SeuratObject::AddMetaData(
    object,
    metadata = as.data.frame(SummarizedExperiment::colData(x), optional = TRUE)
  )
  featureData <- as.data.frame(SummarizedExperiment::rowData(x), optional = TRUE)
  if (ncol(featureData)) {
    object[[assayName]] <- SeuratObject::AddMetaData(
      object[[assayName]],
      metadata = featureData
    )
  }
  coordinates <- as.data.frame(SpatialExperiment::spatialCoords(x))
  coordinates$cell <- colnames(x)
  centroids <- SeuratObject::CreateCentroids(
    coordinates[, c("x", "y", "cell"), drop = FALSE]
  )
  object[[imageName]] <- SeuratObject::CreateFOV(
    coords = list(centroids = centroids),
    type = "centroids",
    molecules = NULL,
    assay = assayName
  )

  if (isTRUE(includeAltExps)) {
    for (alternative in SingleCellExperiment::altExpNames(x)) {
      experiment <- SingleCellExperiment::altExp(x, alternative)
      if (!layer %in% SummarizedExperiment::assayNames(experiment)) {
        warning("Skipping alternative experiment without assay `", layer,
                "`: ", alternative, call. = FALSE)
        next
      }
      alternativeMatrix <- SummarizedExperiment::assay(experiment, layer)
      if (!inherits(alternativeMatrix, "sparseMatrix")) {
        alternativeMatrix <- Matrix::Matrix(alternativeMatrix, sparse = TRUE)
      }
      object[[alternative]] <- SeuratObject::CreateAssay5Object(
        counts = alternativeMatrix
      )
      featureData <- as.data.frame(SummarizedExperiment::rowData(experiment), optional = TRUE)
      if (ncol(featureData)) {
        object[[alternative]] <- SeuratObject::AddMetaData(
          object[[alternative]], metadata = featureData)
      }
    }
  }
  for (name in names(S4Vectors::metadata(x))) {
    value <- S4Vectors::metadata(x)[[name]]
    if (is.data.frame(value)) value <- list(.spamtp_data_frame = value)
    SeuratObject::Misc(object, slot = name) <- value
  }
  labels <- SingleCellExperiment::colLabels(x)
  if (!is.null(labels)) SeuratObject::Idents(object) <- labels
  object
}

#' Convert a Cardinal object to an optional Seurat object
#'
#' @param data A Cardinal `MSImagingExperiment` with aligned spectra.
#' @param multi.run Retained for compatibility; run identities are preserved
#'   during the intermediate conversion.
#' @param seurat.coord Optional replacement coordinates with `X_new` and
#'   `Y_new` columns.
#' @param assay Name of the Seurat assay.
#' @param verbose Show conversion messages.
#'
#' @return A Seurat object.
#' @export
#' @examples
#' utils::str(formals(cardinalToSeurat))
cardinalToSeurat <- function(
    data,
    multi.run = FALSE,
    seurat.coord = NULL,
    assay = "Spatial",
    verbose = TRUE
) {
  if (!is.null(seurat.coord)) {
    pixelData <- Cardinal::pixelData(data)
    pixelData[["x"]] <- seurat.coord$X_new
    pixelData[["y"]] <- seurat.coord$Y_new
    Cardinal::pixelData(data) <- pixelData
  }
  verbose_message(
    "Converting Cardinal through SpatialExperiment to optional Seurat ... ",
    verbose = verbose
  )
  object <- cardinalToSpatialExperiment(data)
  spatialExperimentToSeurat(object, assayName = assay)
}

#' Convert a binned Cardinal matrix to an optional Seurat object
#'
#' @param data Cardinal object supplying coordinates and metadata.
#' @param mtx Binned feature-by-pixel matrix.
#' @param multi.run Retained for compatibility.
#' @param assay Name of the Seurat assay.
#' @param verbose Show conversion messages.
#'
#' @return A Seurat object.
#' @export
#' @examples
#' utils::str(formals(binnedCardinalToSeurat))
binnedCardinalToSeurat <- function(
    data,
    mtx,
    multi.run = FALSE,
    assay = "Spatial",
    verbose = TRUE
) {
  featureData <- Cardinal::MassDataFrame(
    mz = suppressWarnings(as.numeric(rownames(mtx)))
  )
  binned <- Cardinal::MSImagingExperiment(
    spectraData = mtx,
    featureData = featureData,
    pixelData = Cardinal::pixelData(data),
    metadata = S4Vectors::metadata(data)
  )
  cardinalToSeurat(
    binned,
    multi.run = multi.run,
    assay = assay,
    verbose = verbose
  )
}

#' Convert an optional Seurat object to Cardinal
#'
#' @param data A Seurat object.
#' @param assay Seurat assay to convert.
#' @param slot Exact Seurat expression layer name, read with
#'   `SeuratObject::LayerData()`, not an internal S4 slot. This boundary
#'   converter differs from native analysis selectors; see [experimentAccess]
#'   and [seuratToSpatialExperiment()].
#' @param run_col Optional cell-metadata column containing run identifiers.
#' @param feature.metadata Retained for compatibility. Feature metadata are
#'   always preserved.
#' @param verbose Show conversion messages.
#'
#' @return A Cardinal `MSImagingExperiment`.
#' @export
#' @examples
#' utils::str(formals(convertSeuratToCardinal))
convertSeuratToCardinal <- function(
    data,
    assay = "Spatial",
    slot = "counts",
    run_col = NULL,
    feature.metadata = FALSE,
    verbose = TRUE
) {
  verbose_message(
    "Converting optional Seurat input through SpatialExperiment ... ",
    verbose = verbose
  )
  object <- seuratToSpatialExperiment(data, assay = assay, layer = slot)
  if (!is.null(run_col)) {
    columnData <- as.data.frame(SummarizedExperiment::colData(object))
    if (!run_col %in% colnames(columnData)) {
      stop("Run column `", run_col, "` was not found.", call. = FALSE)
    }
    columnData$run <- factor(columnData[[run_col]])
    SummarizedExperiment::colData(object) <- S4Vectors::DataFrame(columnData)
  }
  spatialExperimentToCardinal(object, assayName = slot)
}
