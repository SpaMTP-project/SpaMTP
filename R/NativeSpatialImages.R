.nativeCoordinates <- function(object) {
  object <- .nativeSpatialObject(object)
  coordinates <- SpatialExperiment::spatialCoords(object)
  if (!all(c("x", "y") %in% colnames(coordinates)) ||
      any(!is.finite(coordinates[, c("x", "y"), drop = FALSE])))
    stop("spatialCoords must contain finite x and y coordinates.", call. = FALSE)
  identifiers <- colnames(object)
  if (is.null(identifiers) || anyNA(identifiers) || anyDuplicated(identifiers))
    stop("Spatial pixels must have unique, non-missing names.", call. = FALSE)
  samples <- as.character(SummarizedExperiment::colData(object)$sample_id)
  if (length(samples) != ncol(object) || anyNA(samples))
    stop("colData must contain a sample_id for each pixel.", call. = FALSE)
  data.frame(x = coordinates[, "x"], y = coordinates[, "y"],
             cell = identifiers, sample_id = samples, row.names = identifiers)
}

.nativeImage <- function(object, imageId, sampleId = NULL) {
  images <- SpatialExperiment::imgData(object)
  selected <- which(images$image_id == imageId &
                      (is.null(sampleId) | images$sample_id %in% sampleId))
  if (length(selected) != 1L)
    stop("Select exactly one image using imageId and sampleId.", call. = FALSE)
  row <- images[selected, , drop = FALSE]
  scaleFactor <- row$scaleFactor[[1L]]
  if (!is.finite(scaleFactor) || scaleFactor <= 0)
    stop("The selected image needs a positive imgData scaleFactor.", call. = FALSE)
  raster <- SpatialExperiment::imgRaster(row$data[[1L]])
  list(raster = raster, width = ncol(raster) / scaleFactor,
       height = nrow(raster) / scaleFactor, scaleFactor = scaleFactor,
       sampleId = as.character(row$sample_id[[1L]]))
}

.imageLayer <- function(image, alpha = 1) {
  raster <- image$raster
  if (alpha < 1) {
    colours <- grDevices::adjustcolor(as.vector(raster), alpha.f = alpha)
    # R raster vectors use row-major storage, unlike ordinary matrices.
    raster <- as.raster(matrix(colours, nrow = nrow(raster), byrow = TRUE))
  }
  ggplot2::annotation_raster(raster, xmin = 0, xmax = image$width,
                            ymin = 0, ymax = image$height)
}

.plotNativeSpatialScore <- function(
    object, scores, images = NULL, title = NULL,
    colors = c("darkblue", "lightgrey", "darkred"), guide = "colourbar",
    imageAlpha = 1, pointSize = 1.6, alpha = 1, shape = 16, stroke = 0,
    crop = TRUE, minCutoff = NA, maxCutoff = NA, ncol = NULL,
    interactive = FALSE, imageLabels = NULL
) {
  coordinates <- .nativeCoordinates(object)
  if (length(scores) != nrow(coordinates))
    stop("Scores must contain one value per pixel.", call. = FALSE)
  coordinates$score <- scores
  if (!is.null(images)) {
    entries <- as.data.frame(SpatialExperiment::imgData(object)[,
      c("sample_id", "image_id"), drop = FALSE])
    if (any(!images %in% entries$image_id))
      stop("An image ID was not found in imgData.", call. = FALSE)
    entries <- entries[entries$image_id %in% images, , drop = FALSE]
  } else {
    entries <- data.frame(sample_id = unique(coordinates$sample_id),
                           image_id = NA_character_)
  }
  if (!is.null(imageLabels) && length(imageLabels) != nrow(entries))
    stop("Provide one image label per selected sample/image.", call. = FALSE)
  lower <- .spatialCutoff(scores, minCutoff, -3)
  upper <- .spatialCutoff(scores, maxCutoff, 3)
  plots <- lapply(seq_len(nrow(entries)), function(index) {
    frame <- coordinates[coordinates$sample_id == entries$sample_id[index], ]
    if (!nrow(frame)) stop("The selected image has no pixels.", call. = FALSE)
    label <- title %||% "Pathway score"
    if (nrow(entries) > 1L)
      label <- paste(label, imageLabels[index] %||% entries$sample_id[index])
    plot <- ggplot2::ggplot(frame, ggplot2::aes(x = x, y = y))
    if (!is.na(entries$image_id[index])) {
      image <- .nativeImage(object, entries$image_id[index], entries$sample_id[index])
      plot <- plot + .imageLayer(image, imageAlpha)
      if (!crop) plot <- plot + ggplot2::expand_limits(
        x = c(0, image$width), y = c(0, image$height))
    }
    plot <- plot + ggplot2::geom_point(
      ggplot2::aes(colour = score, fill = score), size = pointSize,
      alpha = alpha, shape = shape, stroke = stroke) +
      ggplot2::scale_colour_gradientn(colours = colors, limits = c(lower, upper),
        oob = scales::squish, guide = guide, name = "z-score") +
      ggplot2::scale_fill_gradientn(colours = colors, limits = c(lower, upper),
        oob = scales::squish, guide = "none") +
      ggplot2::scale_y_reverse() +
      ggplot2::labs(title = label) + ggplot2::theme_void()
    plot <- plot + if (crop) ggplot2::coord_fixed(
      xlim = range(frame$x), ylim = range(frame$y)) else ggplot2::coord_fixed()
    plot
  })
  if (interactive) {
    plots <- lapply(plots, plotly::ggplotly)
    if (length(plots) == 1L) return(plots[[1L]])
    return(plotly::subplot(plots, nrows = ceiling(length(plots) / (ncol %||% 2L))))
  }
  if (length(plots) == 1L) return(plots[[1L]])
  cowplot::plot_grid(plotlist = plots, ncol = ncol)
}
