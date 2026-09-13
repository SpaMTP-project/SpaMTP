.pixelPolygons <- function(coordinates, width, crs = NA_character_) {
  if (is.null(width)) {
    widths <- lapply(split(seq_len(nrow(coordinates)), coordinates$sample_id),
      function(index) {
        if (length(index) < 2L) return(numeric())
        FNN::get.knn(as.matrix(coordinates[index, c("x", "y")]), k = 1)$nn.dist[, 1]
      })
    distances <- unlist(widths)
    width <- stats::median(distances[distances > 0])
  }
  if (length(width) != 1L || !is.finite(width) || width <= 0)
    stop("Supply a positive SM.pixel.width in spatial coordinate units.", call. = FALSE)
  geometry <- lapply(seq_len(nrow(coordinates)), function(index) {
    square <- getSquareCoordinates(coordinates$x[index], coordinates$y[index],
                                    width, coordinates$cell[index])
    sf::st_polygon(list(as.matrix(square[, c("X", "Y")])))
  })
  sf::st_sf(coordinates, geometry = sf::st_sfc(geometry, crs = crs))
}

.mappingWeights <- function(source, target, width, highResolution,
                            radius, polygons, threshold) {
  if (!length(intersect(source$sample_id, target$sample_id)))
    stop("Source and target sample_id values do not match.", call. = FALSE)
  crs <- if (is.null(polygons)) NA_character_ else sf::st_crs(polygons)
  sourcePolygons <- .pixelPolygons(source, width, crs)
  if (highResolution) {
    if (!is.null(polygons))
      stop("ST.polygons is used for area-overlap mapping; set ST.hires = FALSE.",
           call. = FALSE)
    targetGeometry <- sf::st_as_sf(target, coords = c("x", "y"), crs = crs)
  } else if (!is.null(polygons)) {
    if (!inherits(polygons, "sf") || !"cell" %in% colnames(polygons) ||
        anyDuplicated(polygons$cell) || !setequal(polygons$cell, target$cell))
      stop("ST.polygons must be sf polygons with cell IDs matching ST.data.",
           call. = FALSE)
    if (isTRUE(sf::st_is_longlat(polygons)))
      stop("Use planar spatial coordinates and polygons, not longitude/latitude.",
           call. = FALSE)
    targetGeometry <- polygons[match(target$cell, polygons$cell), ]
    if (any(!sf::st_is_valid(targetGeometry)) || any(sf::st_is_empty(targetGeometry)))
      stop("ST.polygons must contain valid, non-empty polygons.", call. = FALSE)
    if (any(!as.character(sf::st_geometry_type(targetGeometry)) %in%
            c("POLYGON", "MULTIPOLYGON")))
      stop("ST.polygons must contain polygon geometries.", call. = FALSE)
  } else {
    if (is.null(radius) || !length(radius) %in% c(1L, nrow(target)) ||
        any(!is.finite(radius)) || any(radius <= 0))
      stop("Supply ST.radius or ST.polygons for area-overlap mapping.", call. = FALSE)
    targetGeometry <- sf::st_buffer(
      sf::st_as_sf(target, coords = c("x", "y"), crs = crs),
      dist = rep(radius, length.out = nrow(target)))
  }
  candidates <- sf::st_intersects(targetGeometry, sourcePolygons)
  matches <- lapply(seq_along(candidates), function(index) {
    selected <- candidates[[index]]
    selected <- selected[source$sample_id[selected] == target$sample_id[index]]
    if (!length(selected)) return(integer())
    if (highResolution) {
      distance <- (source$x[selected] - target$x[index])^2 +
        (source$y[selected] - target$y[index])^2
      return(selected[order(distance, source$cell[selected])][1L])
    }
    targetArea <- as.numeric(sf::st_area(targetGeometry[index, ]))
    if (!is.finite(targetArea) || targetArea <= 0)
      stop("Target polygons must have positive areas.", call. = FALSE)
    coverage <- vapply(selected, function(pixel) {
      intersection <- suppressWarnings(sf::st_intersection(
        sf::st_geometry(targetGeometry[index, ]),
        sf::st_geometry(sourcePolygons[pixel, ])))
      sum(as.numeric(sf::st_area(intersection))) / targetArea
    }, numeric(1))
    selected[coverage > 0 & coverage >= threshold]
  })
  count <- lengths(matches)
  weights <- Matrix::sparseMatrix(
    i = as.integer(unlist(matches, use.names = FALSE)),
    j = rep(seq_along(matches), count),
    x = rep(1 / pmax(count, 1L), count),
    dims = c(nrow(source), nrow(target)), dimnames = list(source$cell, target$cell))
  list(weights = weights, matches = matches)
}

#' Map metabolite intensities to transcriptomic pixels or cells
#'
#' Both inputs must already use the same planar coordinate units and matching
#' sample_id values. Mapping never links different samples. ST.hires uses target
#' centroids inside square MSI pixels; ties use the nearest MSI centre then
#' pixel ID. Otherwise each MSI pixel must cover overlap.threshold of the target
#' spot/polygon area. Matching MSI intensities are averaged with equal weights.
#'
#' @param SM.data,ST.data SpatialExperiment objects with aligned coordinates.
#' @param ST.hires Use point-in-pixel mapping for high-resolution target cells.
#' @param SM.assay,ST.assay Source and target experiment or assay names.
#' @param SM.pixel.width Square MSI pixel width in coordinate units. NULL
#'   estimates the median positive nearest-neighbour distance within samples.
#' @param overlap.threshold Minimum fraction of target area covered by a pixel.
#' @param annotations Preserve MSI rowData and annotation provenance.
#' @param add.metadata Include summaries of source pixel metadata in colData.
#' @param merge.unique.metadata Remove duplicate values in source summaries.
#' @param map.data Map all source assays; FALSE maps counts only.
#' @param new_SPT.assay Name of the transcriptome alternative experiment.
#' @param new_SPM.assay Name of the primary metabolomics experiment.
#' @param verbose Show progress messages.
#' @param ST.radius Target spot radius, scalar or one per target pixel, in
#'   spatial coordinate units; alternatively supply ST.polygons.
#' @param ST.polygons An sf polygon object with a cell column matching ST.data.
#' @param dropUnmapped Remove targets without matching MSI pixels. FALSE keeps
#'   zero-filled intensities and marks these pixels with mapped = FALSE.
#' @return A SpatialExperiment with MSI assays on the target coordinates, target
#'   transcriptomes in altExp, target images and metadata, and mapping provenance.
#' @export
mapSpatialOmics <- function(
    SM.data, ST.data, ST.hires = FALSE, SM.assay = "main", ST.assay = "main",
    SM.pixel.width = NULL, overlap.threshold = 0.2, annotations = TRUE,
    add.metadata = TRUE, merge.unique.metadata = TRUE, map.data = FALSE,
    new_SPT.assay = "transcriptome", new_SPM.assay = "main", verbose = FALSE,
    ST.radius = NULL, ST.polygons = NULL, dropUnmapped = FALSE
) {
  SM.data <- .nativeSpatialObject(SM.data)
  ST.data <- .nativeSpatialObject(ST.data)
  source <- .nativeCoordinates(SM.data)
  target <- .nativeCoordinates(ST.data)
  if (length(overlap.threshold) != 1L || !is.finite(overlap.threshold) ||
      overlap.threshold < 0 || overlap.threshold > 1)
    stop("overlap.threshold must be between zero and one.", call. = FALSE)
  mapping <- .mappingWeights(source, target, SM.pixel.width, ST.hires,
                              ST.radius, ST.polygons, overlap.threshold)
  experiment <- .experimentForAssay(SM.data, SM.assay)
  selected <- if (map.data) SummarizedExperiment::assayNames(experiment) else "counts"
  if (!length(selected) || !all(selected %in% SummarizedExperiment::assayNames(experiment)))
    stop("Source assays are missing; provide counts or set map.data = TRUE.", call. = FALSE)
  assays <- stats::setNames(lapply(selected, function(name) {
    SummarizedExperiment::assay(experiment, name) %*% mapping$weights
  }), selected)
  featureData <- SummarizedExperiment::rowData(experiment)
  if (!annotations) {
    featureData <- S4Vectors::DataFrame(
      mz = .massValues(featureData, rownames(experiment)),
      row.names = rownames(experiment))
  }
  metadata <- SummarizedExperiment::colData(ST.data)
  metadata$SPM_pixels <- S4Vectors::SimpleList(lapply(mapping$matches,
    function(index) source$cell[index]))
  metadata$n_source_pixels <- lengths(mapping$matches)
  metadata$mapped <- lengths(mapping$matches) > 0L
  if (add.metadata) {
    sourceMetadata <- SummarizedExperiment::colData(SM.data)
    for (name in colnames(sourceMetadata)) {
      column <- tail(make.unique(c(colnames(metadata), paste0("SM_", name))), 1L)
      metadata[[column]] <- vapply(mapping$matches, function(index) {
        values <- as.character(sourceMetadata[[name]][index])
        if (merge.unique.metadata) values <- unique(values)
        paste(values, collapse = ", ")
      }, character(1))
    }
  }
  output <- SpatialExperiment::SpatialExperiment(
    assays = assays, rowData = featureData, colData = metadata,
    spatialCoords = SpatialExperiment::spatialCoords(ST.data),
    imgData = SpatialExperiment::imgData(ST.data),
    metadata = S4Vectors::metadata(SM.data))
  SingleCellExperiment::mainExpName(output) <- new_SPM.assay
  SingleCellExperiment::altExp(output, new_SPT.assay) <-
    .experimentForAssay(ST.data, ST.assay)
  S4Vectors::metadata(output)$spamtp_mapping <- list(
    method = if (ST.hires) "centroid within MSI pixel" else "target-area overlap",
    aggregation = "equal-weight mean", overlap_threshold = overlap.threshold,
    source_pixels = source$cell, weights = mapping$weights,
    target_metadata = S4Vectors::metadata(ST.data))
  if (!annotations) {
    S4Vectors::metadata(output)$mz_annotation <- NULL
    S4Vectors::metadata(output)$db_3 <- NULL
  }
  if (dropUnmapped) output <- output[, output$mapped, drop = FALSE]
  verbose_message(paste("Mapped", sum(output$mapped), "of", ncol(output),
                         "target pixels."), verbose = verbose)
  output
}
