# Strict access for migrated workflows: never substitute an unrelated assay.
.nativeSingleSample <- function(object, sampleId = NULL) {
  object <- .nativeSpatialObject(object)
  coordinates <- .nativeCoordinates(object)
  samples <- unique(coordinates$sample_id)
  if (is.null(sampleId)) {
    if (length(samples) != 1L)
      stop("Choose sampleId when multiple samples are present.", call. = FALSE)
    sampleId <- samples
  }
  if (length(sampleId) != 1L || is.na(sampleId) || !sampleId %in% samples)
    stop("Unknown sampleId.", call. = FALSE)
  object[, coordinates$sample_id == sampleId, drop = FALSE]
}

.nativeExpression <- function(object, assay = "main", slot = "counts") {
  experiment <- .experimentForAssay(object, assay)
  if (length(slot) != 1L || is.na(slot) ||
      !slot %in% SummarizedExperiment::assayNames(experiment))
    stop("Expression assay `", slot, "` was not found in `", assay,
         "`. Available: ",
         paste(SummarizedExperiment::assayNames(experiment), collapse = ", "),
         call. = FALSE)
  result <- SummarizedExperiment::assay(experiment, slot)
  if (is.null(rownames(result)) || anyNA(rownames(result)) ||
      anyDuplicated(rownames(result)) || is.null(colnames(result)) ||
      anyNA(colnames(result)) || anyDuplicated(colnames(result)) ||
      !identical(colnames(result), colnames(object)))
    stop("Expression requires unique feature IDs and aligned pixel IDs.",
         call. = FALSE)
  result
}

.mzPlotValues <- function(object, mzs, assay, slot, plusminus = NULL,
                          metabolites = NULL, column = "all_IsomerNames",
                          exact = TRUE) {
  expression <- .nativeExpression(object, assay, slot)
  if (!nrow(expression)) stop("No features to plot.", call. = FALSE)
  features <- rownames(expression)
  masses <- .massValues(.featureMetadata(object, assay), features)
  if (!is.null(plusminus) && (!is.numeric(plusminus) ||
      length(plusminus) != 1L || !is.finite(plusminus) || plusminus < 0))
    stop("plusminus must be one finite non-negative mass tolerance.", call. = FALSE)
  if (is.null(metabolites)) {
    if (!length(mzs) || anyNA(mzs)) stop("Supply m/z values or feature IDs.", call. = FALSE)
    indices <- lapply(mzs, function(target) {
      if (as.character(target) %in% features) return(match(target, features))
      target <- suppressWarnings(as.numeric(target))
      if (!is.finite(target)) stop("Unknown m/z feature.", call. = FALSE)
      which.min(abs(masses - target))
    })
    labels <- paste0("m/z ", format(masses[unlist(indices)], digits = 7, trim = TRUE))
  } else {
    if (!length(metabolites) || anyNA(metabolites))
      stop("Supply metabolite names.", call. = FALSE)
    indices <- lapply(metabolites, function(metabolite) {
      matched <- searchAnnotations(object, metabolite, assay = assay,
                                    column.name = column, search.exact = exact)
      ids <- matched$mz_names %||% rownames(matched)
      index <- match(ids, features)
      if (!length(index) || anyNA(index))
        stop("No matching features for metabolite: ", metabolite, call. = FALSE)
      unique(index)
    })
    labels <- metabolites
  }
  if (!is.null(plusminus)) {
    indices <- lapply(indices, function(centres) {
      unique(unlist(lapply(masses[centres], function(centre) {
        which(abs(masses - centre) <= plusminus +
                8 * .Machine$double.eps * max(1, abs(centre)))
      })))
    })
    labels <- paste0(labels, " +/- ", plusminus)
  }
  values <- do.call(rbind, lapply(indices, function(index) {
    as.numeric(Matrix::colSums(expression[index, , drop = FALSE]))
  }))
  dimnames(values) <- list(make.unique(labels), colnames(expression))
  values
}

.featurePlotEntries <- function(object, coordinates, images, splitBy) {
  if (is.null(images)) {
    entries <- data.frame(sample_id = unique(coordinates$sample_id),
                           image_id = NA_character_)
  } else {
    entries <- as.data.frame(SpatialExperiment::imgData(object)[,
      c("sample_id", "image_id"), drop = FALSE])
    if (!length(images) || any(!images %in% entries$image_id))
      stop("An image ID was not found in imgData.", call. = FALSE)
    entries <- entries[entries$image_id %in% images &
                         entries$sample_id %in% coordinates$sample_id, , drop = FALSE]
  }
  if (!nrow(entries)) stop("No sample/image panels contain selected pixels.", call. = FALSE)
  entries$group <- NA_character_
  if (!is.null(splitBy)) {
    metadata <- .cellMetadata(object)
    if (length(splitBy) != 1L || !splitBy %in% names(metadata))
      stop("split.by must name a colData column.", call. = FALSE)
    coordinates$group <- as.character(metadata[coordinates$cell, splitBy])
    if (anyNA(coordinates$group)) stop("split.by contains missing values.", call. = FALSE)
    entries <- do.call(rbind, lapply(seq_len(nrow(entries)), function(i) {
      groups <- unique(coordinates$group[
        coordinates$sample_id == entries$sample_id[i]])
      entry <- entries[rep(i, length(groups)), , drop = FALSE]
      entry$group <- groups
      entry
    }))
  }
  list(entries = entries, coordinates = coordinates)
}

.combineFeaturePlots <- function(plots, combine, ncol, interactive) {
  if (interactive) plots <- lapply(plots, plotly::ggplotly)
  if (!combine) return(plots)
  if (length(plots) == 1L) return(plots[[1L]])
  if (interactive)
    return(plotly::subplot(plots, nrows = ceiling(length(plots) / (ncol %||% 2L))))
  cowplot::plot_grid(plotlist = plots, ncol = ncol)
}

.plotNativeExpression <- function(
    object, values, images = NULL, sampleIds = NULL, cells = NULL,
    splitBy = NULL, colors = c("lightgrey", "firebrick1"),
    pointSize = 1.6, alpha = 1, imageAlpha = 1, stroke = 0,
    borderColor = "white", pixel = FALSE, crop = TRUE,
    minCutoff = NA, maxCutoff = NA, scale = "feature", axes = FALSE,
    dark = FALSE, coordFixed = TRUE, combine = TRUE, ncol = NULL,
    interactive = FALSE
) {
  coordinates <- .nativeCoordinates(object)
  if (!is.matrix(values) || !nrow(values) ||
      !identical(colnames(values), coordinates$cell) || any(!is.finite(values)))
    stop("Plot values must be finite and aligned by pixel IDs.", call. = FALSE)
  if (!is.null(cells)) {
    if (!length(cells) || any(!cells %in% coordinates$cell))
      stop("Unknown or empty cells selection.", call. = FALSE)
    coordinates <- coordinates[coordinates$cell %in% cells, , drop = FALSE]
  }
  if (!is.null(sampleIds)) {
    if (!length(sampleIds) || any(!sampleIds %in% coordinates$sample_id))
      stop("Unknown or empty sample selection (fov).", call. = FALSE)
    coordinates <- coordinates[coordinates$sample_id %in% sampleIds, , drop = FALSE]
  }
  if (!nrow(coordinates)) stop("No pixels to plot.", call. = FALSE)
  scale <- match.arg(scale, c("feature", "all", "none"))
  if (!length(minCutoff) %in% c(1L, nrow(values)) ||
      !length(maxCutoff) %in% c(1L, nrow(values)))
    stop("Supply one cutoff or one per feature.", call. = FALSE)
  for (cutoff in c(minCutoff, maxCutoff)) {
    if (is.na(cutoff)) next
    quantile <- is.character(cutoff) && grepl("^q[0-9]+([.][0-9]+)?$", cutoff)
    value <- suppressWarnings(as.numeric(if (quantile) sub("^q", "", cutoff) else cutoff))
    if (!is.finite(value) || (quantile && (value < 0 || value > 100)))
      stop("Cutoffs must be finite numbers, NA or q0 through q100.", call. = FALSE)
  }
  if (!length(alpha) %in% c(1L, 2L) || any(!is.finite(alpha)) ||
      any(alpha < 0 | alpha > 1) || length(imageAlpha) != 1L ||
      !is.finite(imageAlpha) || imageAlpha < 0 || imageAlpha > 1)
    stop("Point/image alpha must lie between zero and one.", call. = FALSE)
  selected <- .featurePlotEntries(object, coordinates, images, splitBy)
  coordinates <- selected$coordinates
  entries <- selected$entries
  coordinates <- coordinates[coordinates$sample_id %in% entries$sample_id, , drop = FALSE]
  values <- values[, coordinates$cell, drop = FALSE]
  lower <- upper <- numeric(nrow(values))
  for (i in seq_len(nrow(values))) {
    lower[i] <- .spatialCutoff(values[i, ], rep(minCutoff, length.out = nrow(values))[i],
                              min(values[i, ]))
    upper[i] <- .spatialCutoff(values[i, ], rep(maxCutoff, length.out = nrow(values))[i],
                              max(values[i, ]))
    if (lower[i] > upper[i]) stop("Minimum cutoff exceeds maximum cutoff.", call. = FALSE)
    values[i, ] <- pmax(lower[i], pmin(upper[i], values[i, ]))
  }
  if (scale == "all") {
    lower[] <- min(lower)
    upper[] <- max(upper)
  }
  plots <- list()
  for (i in seq_len(nrow(values))) {
    for (j in seq_len(nrow(entries))) {
      entry <- entries[j, ]
      frame <- coordinates[coordinates$sample_id == entry$sample_id, , drop = FALSE]
      if (!is.na(entry$group)) frame <- frame[frame$group == entry$group, , drop = FALSE]
      frame$score <- as.numeric(values[i, frame$cell])
      limits <- if (scale == "none") range(frame$score) else c(lower[i], upper[i])
      frame$opacity <- rep(alpha[1L], nrow(frame))
      if (length(alpha) == 2L && diff(limits) > 0)
        frame$opacity <- alpha[1L] + diff(alpha) * (frame$score - limits[1L]) / diff(limits)
      labels <- c(rownames(values)[i], if (nrow(entries) > 1L)
        c(entry$sample_id, entry$image_id, entry$group))
      label <- paste(labels[!is.na(labels)], collapse = " | ")
      plot <- ggplot2::ggplot(frame, ggplot2::aes(x = x, y = y))
      if (!is.na(entry$image_id)) {
        image <- .nativeImage(object, entry$image_id, entry$sample_id)
        plot <- plot + .imageLayer(image, imageAlpha)
        if (!crop) plot <- plot + ggplot2::expand_limits(
          x = c(0, image$width), y = c(0, image$height))
      }
      plot <- plot + ggplot2::geom_point(
        ggplot2::aes(fill = score, alpha = opacity),
        shape = if (pixel) 22 else 21, colour = borderColor,
        size = pointSize, stroke = stroke) +
        ggplot2::scale_alpha_identity() +
        ggplot2::scale_fill_gradientn(colours = colors, limits = limits,
          oob = scales::squish, name = rownames(values)[i]) +
        ggplot2::scale_y_reverse() + ggplot2::labs(title = label)
      plot <- plot + if (axes) ggplot2::theme_minimal() else ggplot2::theme_void()
      if (dark) plot <- plot + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "black", colour = NA),
        text = ggplot2::element_text(colour = "white"))
      coord <- if (coordFixed) ggplot2::coord_fixed else ggplot2::coord_cartesian
      plot <- plot + if (crop) coord(xlim = range(frame$x), ylim = range(frame$y)) else coord()
      plots[[length(plots) + 1L]] <- plot
    }
  }
  .combineFeaturePlots(plots, combine, ncol, interactive)
}
