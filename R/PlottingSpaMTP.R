# MSI plotting and feature queries for native Bioconductor containers.

#' Finds the nearest m/z peak to a given value in a SpaMTP Object
#'
#' @param data A SummarizedExperiment, including SingleCellExperiment or
#'   SpatialExperiment, with feature m/z values in `rowData()` or feature names.
#' @param target_mz Numeric value defining the target m/z peak
#' @param assay Primary (`"main"` or NULL) or alternative experiment name.
#'   The mass lookup uses the selected experiment's `rowData()`, not an
#'   expression layer; see [experimentAccess].
#'
#' @returns The exact character feature identifier for the nearest measured
#'   m/z, not a numeric mass. Ties select the first feature in the stored order.
#' @export
#'
#' @examples
#' utils::str(formals(findNearestMZ))
#' # findNearestMZ(spe, target_mz = 400.01)
findNearestMZ <- function(data, target_mz, assay = NULL){
  featureMetadata <- .featureMetadata(data, assay)
  featureNames <- rownames(featureMetadata)
  if (is.null(featureNames)) {
    featureNames <- rownames(.assayData(data, assay, "counts"))
  }
  numbers <- .massValues(featureMetadata, featureNames)
  if (length(target_mz) != 1L || !is.finite(target_mz)) {
    stop("`target_mz` must be one finite numeric value.", call. = FALSE)
  }
  featureNames[[which.min(abs(numbers - target_mz))]]
}


#' Sum selected metabolite features per pixel
#'
#' The sum is stored in colData without modifying the expression assays.
#' @param data A SpatialExperiment or aligned Cardinal experiment.
#' @param mzs Character feature IDs to sum. Duplicates count once.
#' @param assay Primary or alternative experiment.
#' @param slot Expression assay.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param bin_name Output colData column name. Existing columns are protected.
#' @return The updated SpatialExperiment.
#' @export
#' @examples
#' utils::str(formals(binMetabolites))
binMetabolites <- function(data, mzs, assay = "main", slot = "counts",
                           bin_name = "Binned_Metabolites") {
  data <- .nativeSpatialObject(data)
  values <- .nativeExpression(data, assay, slot)
  if (!length(mzs) || anyNA(mzs) || any(!mzs %in% rownames(values)))
    stop("Supply existing feature IDs.", call. = FALSE)
  if (length(bin_name) != 1L || is.na(bin_name) || !nzchar(bin_name) ||
      bin_name %in% colnames(SummarizedExperiment::colData(data)))
    stop("Choose a new colData column name.", call. = FALSE)
  SummarizedExperiment::colData(data)[[bin_name]] <-
    as.numeric(Matrix::colSums(values[unique(mzs), , drop = FALSE]))
  data
}

#' Convert point symbols to squares without changing coordinates or scales
#'
#' For combined feature panels, request combine=FALSE from the original plot
#' function and pass the resulting list. This function changes point symbols,
#' not tissue geometry or the physical pixel footprint.
#' @param plot A ggplot with point layers, or a list of such plots.
#' @return The ggplot or list, retaining images, coordinates, scales and labels.
#' @export
#' @examples
#' utils::str(formals(pixelPlot))
pixelPlot <- function(plot) {
  if (is.list(plot) && !inherits(plot, "ggplot")) return(lapply(plot, pixelPlot))
  if (!inherits(plot, "ggplot")) stop("Supply a ggplot or list of ggplots.", call. = FALSE)
  points <- which(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomPoint"),
                          logical(1)))
  if (!length(points))
    stop("No point layers; use combine=FALSE for multi-panel feature plots.", call. = FALSE)
  for (i in points) {
    layer <- plot$layers[[i]]
    filled <- !is.null(layer$mapping$fill %||% plot$mapping$fill)
    parameters <- utils::modifyList(layer$aes_params, list(shape = if (filled) 22 else 15))
    plot$layers[[i]] <- ggplot2::ggproto(NULL, layer, aes_params = parameters)
  }
  plot
}




#' Plot MSI features from a Bioconductor spatial container
#'
#' These functions read the primary experiment or a named altExp without
#' modifying the object. Numeric m/z queries select the nearest measured mass;
#' character feature IDs select exact rows. A mass tolerance sums all peaks
#' within that distance of the selected measured mass. Overlapping annotation
#' windows count each peak once. Multiple samples are plotted separately.
#'
#' @param object A SpatialExperiment or aligned Cardinal MSImagingExperiment.
#' @param mzs Numeric m/z values or character feature IDs.
#' @param plusminus Optional non-negative absolute m/z tolerance.
#' @param fov Optional sample IDs from colData(object)$sample_id.
#' @param cols Gradient colours.
#' @param size Point size.
#' @param min.cutoff,max.cutoff One cutoff or one per query; numeric or
#'   quantile strings such as "q10". NA uses the observed range.
#' @param split.by Optional colData column for separate panels within samples.
#' @param alpha Point opacity, or a two-value range mapped to expression.
#' @param border.color,border.size Point border colour and width.
#' @param dark.background Use a black plot background.
#' @param crop Crop to the selected pixels.
#' @param cells Optional pixel IDs to display.
#' @param scale Colour limits: "feature" across samples, "all" across all
#'   queries, or "none" separately for each panel.
#' @param axes Show axes.
#' @param combine Combine panels; FALSE returns a list of ggplots.
#' @param coord.fixed Use a fixed coordinate aspect ratio.
#' @param assay Primary experiment ("main") or altExp name.
#' @param slot Expression assay name, e.g. "counts" or "logcounts".
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param plot.pixel Draw square rather than circular point symbols.
#' @param verbose Retained for call compatibility.
#' @param ... Former Seurat segmentation, molecule and blending arguments
#'   are unsupported and cause an error; convert and plot these separately.
#' @return A ggplot, or a list of panels when combine is FALSE.
#' @export
#' @examples
#' x <- SpatialExperiment::SpatialExperiment(
#'   assays = list(counts = rbind(a = 1:4, b = 4:1)),
#'   rowData = S4Vectors::DataFrame(mz = c(100, 200)),
#'   colData = S4Vectors::DataFrame(row.names = paste0("p", 1:4)),
#'   spatialCoords = cbind(x = 1:4, y = c(0, 1, 0, 1)))
#' imageMZPlot(x, mzs = 100)
imageMZPlot <- function(
    object, mzs, plusminus = NULL, fov = NULL,
    cols = c("lightgrey", "firebrick1"), size = 0.5,
    min.cutoff = NA, max.cutoff = NA, split.by = NULL, alpha = 1,
    border.color = "white", border.size = 0, dark.background = TRUE,
    crop = FALSE, cells = NULL, scale = c("feature", "all", "none"),
    axes = FALSE, combine = TRUE, coord.fixed = TRUE,
    assay = "main", slot = "counts", plot.pixel = FALSE, verbose = TRUE, ...
) {
  if (length(list(...))) stop("Seurat-specific plotting arguments are unsupported.", call. = FALSE)
  object <- .nativeSpatialObject(object)
  values <- .mzPlotValues(object, mzs, assay, slot, plusminus)
  .plotNativeExpression(object, values, sampleIds = fov, colors = cols,
    pointSize = size, minCutoff = min.cutoff, maxCutoff = max.cutoff,
    splitBy = split.by, alpha = alpha, borderColor = border.color,
    stroke = border.size %||% 0, dark = dark.background, crop = crop,
    cells = cells, scale = match.arg(scale), axes = axes, combine = combine,
    coordFixed = coord.fixed, pixel = plot.pixel)
}

#' Plot annotated metabolites without an optical image
#'
#' All matching peaks are summed per metabolite. A peak shared by overlapping
#' mass-tolerance windows is counted once. No temporary metadata is written.
#' @inheritParams imageMZPlot
#' @param metabolites Metabolite names to search in feature annotations.
#' @param column.name Annotation column in rowData.
#' @param plot.exact Match complete semicolon-delimited names (case-insensitive).
#' @return A ggplot or list of panels.
#' @export
#' @examples
#' utils::str(formals(imageMZAnnotationPlot))
imageMZAnnotationPlot <- function(
    object, metabolites, plusminus = NULL, fov = NULL,
    cols = c("lightgrey", "firebrick1"), size = 0.5,
    min.cutoff = NA, max.cutoff = NA, split.by = NULL, alpha = 1,
    border.color = "white", border.size = 0, dark.background = TRUE,
    crop = FALSE, cells = NULL, scale = c("feature", "all", "none"),
    axes = FALSE, combine = TRUE, coord.fixed = TRUE,
    assay = "main", slot = "counts", column.name = "all_IsomerNames",
    plot.exact = TRUE, plot.pixel = FALSE, verbose = TRUE, ...
) {
  if (length(list(...))) stop("Seurat-specific plotting arguments are unsupported.", call. = FALSE)
  object <- .nativeSpatialObject(object)
  values <- .mzPlotValues(object, NULL, assay, slot, plusminus,
    metabolites, column.name, plot.exact)
  .plotNativeExpression(object, values, sampleIds = fov, colors = cols,
    pointSize = size, minCutoff = min.cutoff, maxCutoff = max.cutoff,
    splitBy = split.by, alpha = alpha, borderColor = border.color,
    stroke = border.size %||% 0, dark = dark.background, crop = crop,
    cells = cells, scale = match.arg(scale), axes = axes, combine = combine,
    coordFixed = coord.fixed, pixel = plot.pixel)
}

#' Plot MSI features with optional optical images
#'
#' Images come from imgData and use its scaleFactor to place the raster in
#' spatialCoords units. Each sample/image is a separate panel. With images=NULL
#' all samples are plotted without an image.
#' @inheritParams imageMZPlot
#' @param images Optional image IDs in imgData.
#' @param keep.scale Colour scaling: "feature", "all" or "none"; NULL means "none".
#' @param ncol Number of columns in combined plots.
#' @param pt.size.factor Point size.
#' @param image.alpha Optical image opacity.
#' @param stroke Point border width.
#' @param interactive Return Plotly widgets rather than static plots.
#' @param information Deprecated Seurat hover data; non-NULL is unsupported.
#' @return A ggplot or Plotly widget; a list when combine=FALSE.
#' @export
#' @examples
#' utils::str(formals(spatialMZPlot))
spatialMZPlot <- function(
    object, mzs, plusminus = NULL, images = NULL, crop = TRUE,
    assay = "main", slot = "counts", keep.scale = "feature",
    min.cutoff = NA, max.cutoff = NA, ncol = NULL, combine = TRUE,
    pt.size.factor = 1.6, alpha = c(1, 1), image.alpha = 1,
    stroke = 0.25, interactive = FALSE, information = NULL, verbose = TRUE
) {
  if (!is.null(information))
    stop("information is a retired Seurat hover-data argument.", call. = FALSE)
  object <- .nativeSpatialObject(object)
  values <- .mzPlotValues(object, mzs, assay, slot, plusminus)
  .plotNativeExpression(object, values, images = images, crop = crop,
    scale = keep.scale %||% "none", minCutoff = min.cutoff,
    maxCutoff = max.cutoff, ncol = ncol, combine = combine,
    pointSize = pt.size.factor, alpha = alpha, imageAlpha = image.alpha,
    stroke = stroke, interactive = interactive)
}

#' Plot annotated metabolites with optional optical images
#'
#' Matching peaks are summed using the same rules as imageMZAnnotationPlot.
#' @inheritParams spatialMZPlot
#' @inheritParams imageMZAnnotationPlot
#' @return A ggplot or Plotly widget; a list when combine=FALSE.
#' @export
#' @examples
#' utils::str(formals(spatialMZAnnotationPlot))
spatialMZAnnotationPlot <- function(
    object, metabolites, plusminus = NULL, images = NULL, crop = TRUE,
    assay = "main", slot = "counts", keep.scale = "feature",
    min.cutoff = NA, max.cutoff = NA, ncol = NULL, combine = TRUE,
    pt.size.factor = 1.6, alpha = c(1, 1), image.alpha = 1,
    stroke = 0.25, interactive = FALSE, information = NULL,
    column.name = "all_IsomerNames", plot.exact = TRUE, verbose = TRUE
) {
  if (!is.null(information))
    stop("information is a retired Seurat hover-data argument.", call. = FALSE)
  object <- .nativeSpatialObject(object)
  values <- .mzPlotValues(object, NULL, assay, slot, plusminus,
    metabolites, column.name, plot.exact)
  .plotNativeExpression(object, values, images = images, crop = crop,
    scale = keep.scale %||% "none", minCutoff = min.cutoff,
    maxCutoff = max.cutoff, ncol = ncol, combine = combine,
    pointSize = pt.size.factor, alpha = alpha, imageAlpha = image.alpha,
    stroke = stroke, interactive = interactive)
}


#' Plot mean mass spectra from a native spatial experiment
#'
#' Feature masses come from rowData, so arbitrary feature IDs are supported.
#' @param data A SpatialExperiment or aligned Cardinal experiment.
#' @param group.by Optional colData column defining overlaid groups.
#' @param split.by Optional colData column defining separate panels; mutually
#'   exclusive with group.by.
#' @param cols Optional group colours.
#' @param assay Primary or alternative experiment.
#' @param slot Expression assay.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param label.annotations Display annotations instead of m/z labels.
#' @param annotation.column Annotation column in rowData.
#' @param mz.labels Numeric m/z values to label using nearest measured peaks.
#' @param metabolite.labels Metabolite names to label; mutually exclusive with mz.labels.
#' @param xlab,ylab Axis labels.
#' @param mass.range Optional two-value displayed mass range.
#' @param y.lim Optional y-axis limits.
#' @param labelCex Text size.
#' @param labelAdj,labelOffset Vertical and horizontal text justification.
#' @param labelCol Text colour.
#' @param nlabels.to.show Maximum annotations displayed per feature.
#' @return A ggplot.
#' @export
#' @examples
#' utils::str(formals(massIntensityPlot))
massIntensityPlot <- function(
    data, group.by = NULL, split.by = NULL, cols = NULL, assay = "main",
    slot = "counts", label.annotations = FALSE, annotation.column = "all_IsomerNames",
    mz.labels = NULL, metabolite.labels = NULL, xlab = "m/z", ylab = "intensity",
    mass.range = NULL, y.lim = NULL, labelCex = 5, labelAdj = -1,
    labelOffset = 0, labelCol = "#eb4034", nlabels.to.show = NULL
) {
  data <- .nativeSpatialObject(data)
  expression <- .nativeExpression(data, assay, slot)
  metadata <- .featureMetadata(data, assay)
  masses <- .massValues(metadata, rownames(expression))
  if (!is.null(group.by) && !is.null(split.by))
    stop("Choose group.by or split.by, not both.", call. = FALSE)
  if (!is.null(mz.labels) && !is.null(metabolite.labels))
    stop("Choose mz.labels or metabolite.labels, not both.", call. = FALSE)
  grouping <- group.by %||% split.by
  if (is.null(grouping)) {
    groups <- rep("Sample mean", ncol(data))
  } else {
    if (length(grouping) != 1L || !grouping %in% colnames(SummarizedExperiment::colData(data)))
      stop("Grouping column was not found in colData.", call. = FALSE)
    groups <- as.character(SummarizedExperiment::colData(data)[[grouping]])
    if (anyNA(groups)) stop("Grouping contains missing values.", call. = FALSE)
  }
  selected <- integer()
  if (!is.null(mz.labels)) {
    if (!is.numeric(mz.labels) || any(!is.finite(mz.labels)))
      stop("mz.labels must be finite numeric masses.", call. = FALSE)
    selected <- vapply(mz.labels, function(mass) which.min(abs(masses - mass)), integer(1))
  }
  if (!is.null(metabolite.labels)) {
    selected <- unique(unlist(lapply(metabolite.labels, function(metabolite) {
      rows <- searchAnnotations(data, metabolite, assay, search.exact = TRUE,
                                column.name = annotation.column)
      if (!nrow(rows)) stop("No annotation for: ", metabolite, call. = FALSE)
      match(rows$mz_names, rownames(expression))
    })))
  }
  labels <- rep(NA_character_, nrow(expression))
  if (length(selected)) {
    if (label.annotations) {
      if (is.null(annotation.column) || !annotation.column %in% names(metadata))
        stop("Annotation column was not found in rowData.", call. = FALSE)
      labels[selected] <- as.character(metadata[[annotation.column]][selected])
      if (!is.null(nlabels.to.show))
        labels[selected] <- labels_to_show(labels[selected], n = nlabels.to.show)
    } else labels[selected] <- paste0("m/z ", format(masses[selected], trim = TRUE))
  }
  frames <- lapply(unique(groups), function(group) {
    data.frame(mass = masses, intensity = as.numeric(Matrix::rowMeans(
      expression[, groups == group, drop = FALSE])),
      variable = group, annotation = labels)
  })
  frame <- do.call(rbind, frames)
  if (any(!is.finite(frame$intensity))) stop("Expression must be finite.", call. = FALSE)
  if (!is.null(mass.range)) {
    if (length(mass.range) != 2L || any(!is.finite(mass.range)) || diff(mass.range) < 0)
      stop("mass.range must contain increasing finite limits.", call. = FALSE)
    frame <- frame[frame$mass >= mass.range[1L] & frame$mass <= mass.range[2L], ]
  }
  annotations <- frame[!is.na(frame$annotation), , drop = FALSE]
  if (is.null(split.by)) {
    annotations <- annotations[order(annotations$mass, -annotations$intensity), ]
    annotations <- annotations[!duplicated(annotations$mass), ]
  }
  plot <- ggplot2::ggplot(frame, ggplot2::aes(x = mass, y = intensity, colour = variable)) +
    ggplot2::geom_line() + ggplot2::theme_classic() +
    ggplot2::labs(x = xlab, y = ylab, colour = grouping %||% "Group")
  if (nrow(annotations)) plot <- plot + ggplot2::geom_text(
    data = annotations, ggplot2::aes(label = annotation), size = labelCex,
    vjust = labelAdj, hjust = labelOffset, colour = labelCol)
  if (!is.null(split.by)) plot <- plot + ggplot2::facet_wrap(~variable, ncol = 1, scales = "free_y")
  if (!is.null(cols)) plot <- plot + ggplot2::scale_colour_manual(values = cols)
  if (!is.null(y.lim)) plot <- plot + ggplot2::coord_cartesian(ylim = y.lim)
  plot
}



#' Compare two spatial features in 3D
#'
#' Expression comes from the primary or alternative experiment; numeric
#' colData columns can also be plotted. Coordinates and images use
#' SpatialExperiment accessors. Only one sample is displayed at a time.
#' @param data A SpatialExperiment or aligned Cardinal experiment.
#' @param features One or two feature IDs or numeric colData column names.
#'   A single feature is repeated for comparison across two experiments.
#' @param assays One or two primary/altExp names.
#' @param slots One or two expression assay names.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param between.layer.height Distance between layers.
#' @param names Optional layer labels.
#' @param size Marker size.
#' @param col.palette One or two Plotly colour palettes.
#' @param x.axis.label,y.axis.label,z.axis.label Axis labels.
#' @param show.x.ticks,show.y.ticks,show.z.ticks Show axis tick labels.
#' @param show.image Optional image ID from imgData.
#' @param plot.height,plot.width Widget dimensions.
#' @param image.sf Retired Seurat scale-factor selector; must be NULL.
#'   The selected imgData row supplies the scaleFactor.
#' @param downscale.image Optional positive integer raster sampling stride.
#' @param sampleId Sample ID; required when data contains multiple samples.
#' @return A Plotly widget.
#' @export
#' @rawNamespace import(plotly, except = last_plot)
#' @examples
#' utils::str(formals(plot3DFeature))
plot3DFeature <- function(
    data, features, assays = c("transcriptome", "main"), slots = "counts",
    between.layer.height = 100, names = NULL, size = 3, col.palette = "Reds",
    x.axis.label = "x", y.axis.label = "y", z.axis.label = "z",
    show.x.ticks = FALSE, show.y.ticks = FALSE, show.z.ticks = FALSE,
    show.image = NULL, plot.height = 800, plot.width = 1500,
    image.sf = NULL, downscale.image = NULL, sampleId = NULL
) {
  data <- .nativeSingleSample(data, sampleId)
  if (!is.null(image.sf)) stop("Use imgData scaleFactor instead of image.sf.", call. = FALSE)
  arguments <- list(features, assays, slots, col.palette)
  if (any(!lengths(arguments) %in% c(1L, 2L)))
    stop("Supply one or two features, assays, slots and palettes.", call. = FALSE)
  features <- rep(features, length.out = 2L)
  assays <- rep(assays, length.out = 2L)
  slots <- rep(slots, length.out = 2L)
  palettes <- rep(as.list(col.palette), length.out = 2L)
  if (!is.null(names) && !length(names) %in% c(1L, 2L))
    stop("Supply one or two layer labels.", call. = FALSE)
  labels <- if (is.null(names)) paste(features, assays, sep = " | ") else
    rep(names, length.out = 2L)
  coordinates <- .nativeCoordinates(data)
  metadata <- .cellMetadata(data)
  plot <- plotly::plot_ly(height = plot.height, width = plot.width)
  for (i in seq_len(2L)) {
    if (features[i] %in% colnames(metadata)) {
      values <- metadata[[features[i]]]
    } else {
      expression <- .nativeExpression(data, assays[i], slots[i])
      if (!features[i] %in% rownames(expression)) stop("Unknown feature: ", features[i], call. = FALSE)
      values <- as.numeric(expression[features[i], ])
    }
    if (!is.numeric(values) || any(!is.finite(values)))
      stop("Features must contain finite numeric values.", call. = FALSE)
    plot <- plotly::add_trace(plot, x = coordinates$x, y = coordinates$y,
      z = rep((2L - i) * between.layer.height, nrow(coordinates)),
      text = coordinates$cell, type = "scatter3d", mode = "markers",
      name = labels[i], marker = list(color = values, size = size,
        colorscale = palettes[[i]], showscale = TRUE,
        colorbar = list(x = (i - 1L) * 1.1, title = list(text = labels[i]))))
  }
  if (!is.null(show.image)) {
    image <- .nativeImage(data, show.image, unique(data$sample_id))
    stride <- downscale.image %||% 1L
    if (length(stride) != 1L || !is.finite(stride) || stride < 1 || stride != floor(stride))
      stop("downscale.image must be a positive integer.", call. = FALSE)
    rows <- seq.int(1L, nrow(image$raster), by = stride)
    columns <- seq.int(1L, ncol(image$raster), by = stride)
    # Match the row-major colour vector returned by the raster subset.
    grid <- expand.grid(column = columns, row = rows)
    colours <- as.vector(image$raster[rows, columns, drop = FALSE])
    plot <- plotly::add_trace(plot,
      x = (grid$column - 0.5) / image$scaleFactor,
      y = (grid$row - 0.5) / image$scaleFactor,
      z = rep(-between.layer.height, nrow(grid)),
      type = "scatter3d", mode = "markers", name = show.image,
      marker = list(color = colours, size = size, showscale = FALSE))
  }
  plotly::layout(plot, scene = list(aspectmode = "data",
    xaxis = list(title = x.axis.label, showticklabels = show.x.ticks),
    yaxis = list(title = y.axis.label, showticklabels = show.y.ticks,
                 autorange = "reversed"),
    zaxis = list(title = z.axis.label, showticklabels = show.z.ticks)))
}




#' Export an interactive spatial detection-density viewer
#'
#' The spatial kernel uses coordinates of nonzero pixels, not intensity
#' weights. The mass-spectrum overlay uses an intensity-weighted density
#' without replicating each mass according to its intensity. This exploratory
#' viewer is not a spatial significance test. Its HTML requires internet access
#' to load the existing Plotly and JavaScript KDE libraries.
#' @param object A SpatialExperiment or aligned Cardinal experiment.
#' @param assay Primary or alternative experiment.
#' @param slot Non-negative intensity assay.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param folder Output directory.
#' @param sampleId Sample ID; required for multiple samples.
#' @param ... Reserved; unused arguments cause an error.
#' @return Invisibly, the path to mzs_density_map.html.
#' @export
#' @examples
#' utils::str(formals(densityMap))
densityMap <- function(object, assay = "main", slot = "counts",
                        folder = getwd(), sampleId = NULL, ...) {
  if (length(list(...))) stop("Unused densityMap arguments.", call. = FALSE)
  object <- .nativeSingleSample(object, sampleId)
  expression <- .nativeExpression(object, assay, slot)
  if (!nrow(expression)) stop("No features to export.", call. = FALSE)
  annotated_table <- .featureMetadata(object, assay)
  masses <- .massValues(annotated_table, rownames(expression))
  indices <- .nativeCoordinates(object)[, c("x", "y"), drop = FALSE]
  means <- as.numeric(Matrix::rowMeans(expression))
  if (any(!is.finite(means)) || any(means < 0))
    stop("Density export requires finite non-negative intensities.", call. = FALSE)
  labels <- annotated_table$all_IsomerNames %||%
    annotated_table$annotation %||% rep("", nrow(expression))
  identifiers <- annotated_table$all_Isomers_IDs %||% rep("", nrow(expression))
  annotationRows <- lapply(seq_len(nrow(expression)), function(i) {
    list(masses[i], rownames(expression)[i], labels[i], identifiers[i])
  })
  # Escape HTML delimiters as well as JSON strings before embedding in script.
  scriptJSON <- function(value, ...) {
    json <- jsonlite::toJSON(value, digits = NA, na = "null", ...)
    gsub("<", "\\u003c", json, fixed = TRUE)
  }
  annotated_json <- scriptJSON(annotationRows)
  histogram_data_to_be_added <- scriptJSON(
    data.frame(mz = masses, mean = means), dataframe = "rows")
  if (sum(means) > 0) {
    kde <- stats::density(masses, weights = means / sum(means), bw = 0.05)
    kde$y <- kde$y * max(means) / max(kde$y)
  } else {
    kde <- list(x = seq(min(masses) - 0.15, max(masses) + 0.15, length.out = 512),
                y = rep(0, 512))
  }
  kde_json <- scriptJSON(list(kde$x, kde$y))
  locations <- lapply(seq_len(nrow(expression)), function(i) {
    values <- as.numeric(expression[i, ])
    if (any(!is.finite(values)) || any(values < 0))
      stop("Density export requires finite non-negative intensities.", call. = FALSE)
    unname(as.matrix(indices[values > 0, , drop = FALSE]))
  })
  json_array <- scriptJSON(locations)


  javascript <- paste0(
    "
<!DOCTYPE html>
  <div id='Data'></div>
    <html lang='en'>

      <head>
      <meta charset='UTF-8'>
        <meta name='viewport' content='width=device-width, initial-scale=1.0'>
          <title>Interactive Bar and Dot Plot</title>
          <script src='https://cdn.plot.ly/plotly-latest.min.js'></script>
            <script src='https://cdnjs.cloudflare.com/ajax/libs/PapaParse/4.1.2/papaparse.min.js'></script>
              <script src='https://d3js.org/d3.v7.min.js'></script>
                  <script type='text/javascript' src='https://cdn.jsdelivr.net/gh/stdlib-js/stats-kde2d@umd/browser.js'></script>
                    <script type='text/javascript' src='https://cdn.jsdelivr.net/gh/stdlib-js/ndarray@umd/browser.js'></script>
                      </head>

                      <body>
                      <script type='text/javascript'>
                        (function () {
                          window.kde2d;
                        })();
                      </script>
                        <script type='text/javascript'>
                          (function () {
                            window.ndarray;
                          })();
                        </script>

    <style>
      #container {
          display: flex;
          flex-direction: column;
      }
      #plots {
          display: flex;
      }
      #barPlot, #secPlot {
          width: 100%;
          border: 0px solid black; /* Optional: to visualize the divs */
          box-sizing: border-box; /* Ensure borders are included in the dimensions */
          height: 450px;
          padding: 10px;
          margin-bottom: 0px;
          z-index: 1;
      }
      #buttonContainer {
        display: flex;
        flex-direction: column;
        position: absolute;
        left: 55%; /* Adjust to move horizontally within secPlot */
        top: 5px; /* Adjust to move vertically within secPlot */
      }
      #tablePlot {
          width: 100%;
          border: 2px solid black; /* Optional: to visualize the div */
          box-sizing: border-box; /* Ensure borders are included in the dimensions */
          height: 300px;
          padding: 30px;
      }
      #message {
        display: flex;
        flex-direction: column;
        position: absolute;
        left: 55%; /* Adjust to move horizontally within secPlot */
        top: 50px; /* Adjust to move vertically within secPlot */
      }
  </style>
      <div id='container'>
        <div id='plots'>
            <div id='barPlot'></div>
            <div id='secPlot'></div>
              <div id = 'warning'></div>
              <div id='buttonContainer'>
              <div id='message'>
                  <strong>Instructions:</strong> Entry a number to describe the number of equally spaced points at which the density is to be estimated.
                  <input type='number' id='numberEntry' min='2' max='100' width='30' value='30'/>
              </div>
              </div>
        </div>
        <div id='tablePlot'>
        </div>
    </div>
  <script>
    const warning = document.getElementById('warning');
    document.getElementById('numberEntry').addEventListener('input', function() {
    var value = parseInt(this.value);
    if (value > 100) {
        this.value = 100; // Cap the value at 100
    }
});
          const data = ",
    histogram_data_to_be_added,
    ";
          // Add event listener for unhover on bar plot
          var json_array =",
    json_array,
    "
          var annotated_table =",
    annotated_json,
    "
          var kde =",
    kde_json,
    "
          var trace = {
              x: kde[0],
              y: kde[1],
              mode: 'lines',
              name: 'Estimated average density spectrum',
              line: {shape: 'spline'},
              type: 'scatter',
              width: 0.3,
              hoverinfo: null
              //hoverinfo: 'Estimated density' + kde[1]
            };

            const xValues = data.map(d => d.mz);
            const yValues = data.map(d => d.mean);

            // Create trace for the bar plot
            var barTrace = {
              x: xValues,
              y: yValues,
              type: 'bar',
              width: 0.3,
              name: 'Peak bins',
              hoverinfo: 'm/z value' + xValues
            };

            // Create trace for the dot plot


            // Set layout for bar plot
            var barLayout = {
              hovermode: 'x unified',
              title: 'Bar Plot',
              xaxis: {
                title: 'm/z values'
              },
              yaxis: {
                title: 'intensity'
              },
              plot_bgcolor: 'rgba(0,0,0,0)', // Transparent background for the plot area
              paper_bgcolor: 'rgba(0,0,0,0)',
              hoverdistance: 50
            };

            // Set layout for second plot

            // Create bar plot
            Plotly.newPlot('barPlot', [trace,barTrace], barLayout);

            function findArrayByNumber(array, number) {
              let searchString = 'mz-' + number;
              for (let i = 0; i < array.length; i++) {
                let element = array[i][0];
                if (element === searchString) {
                  return i;
                }
              }
              return -1;
            }
            document.querySelectorAll('.plotly .scatterlayer .trace').forEach(function(el, index) {
            if(index === 0) {
            el.style.pointerEvents = 'none';  // Disable interactivity for this trace
            }
            });
            // Add event listener for hover on bar plot
            var ind = 0
            function hideMessage() {
              warning.textContent = ''; // Clear the message text
              warning.style.display = 'none'; // Hide the message box
        }
            document.getElementById('barPlot').on('plotly_click', function (data_event){
              hideMessage()
              if(data_event.points.length == 2){
                var pointNumber = data_event.points[1].pointNumber;
              }else if (data_event.points[0].fullData.type == 'bar'){
                var pointNumber = data_event.points[0].pointNumber;
              }else{
                Plotly.purge('secPlot');
                warning.textContent = 'No explicit peak bin selected, consider closing the density spetrum and zoom in to select explicit m/z peak'; // Set the message text
                warning.style.display = 'block'; // Make the message box visible
                return;
              }

              ind = pointNumber
              console.log(pointNumber)
              var selected_mz = xValues[pointNumber]
              Plotly.purge('secPlot');
              displaydensity(pointNumber)
              var table_index = pointNumber;
              if(table_index != -1){
                displayTable(table_index);
              }else{
                var table_update = [[selected_mz],['Not annotated by given adduct/db'],
                                    ['Not annotated by given adduct/db'],
                                    ['Not annotated by given adduct/db']]
                var tableData = [
                  {
                    type: 'table',
                    header: {
                      values: ['mz', 'mz_name', 'annotated metabolites', 'entry of the metabolites'],
                      align: 'center',
                      line: {width: 1, color: 'black'},
                      fill: {color: 'grey'},
                      font: {family: 'Arial', size: 12, color: 'white'}
                    },
                    cells: {
                      values: table_update,
                      align: 'center',
                      line: {color: 'black', width: 1},
                      fill: {color: ['white', 'lightgrey']},
                      font: {family: 'Arial', size: 11, color: ['black']}
                    }
                  }
                ];
                var layout_table = {
                  title: 'Annotation table for m/z: ' + data[table_index].mz,
                  autosize: true,
                  plot_bgcolor: 'rgba(0,0,0,0)', // Transparent background for the plot area
                  paper_bgcolor: 'rgba(0,0,0,0)',
                };
                Plotly.newPlot('tablePlot', tableData,layout_table);
              }
            });
            // Get the density plot done
            // If the element exists, hide the loading message
            var space  = 30
            document.getElementById('numberEntry').addEventListener('blur', function(event) {
            space = Math.max(2, Math.min(100, Number(event.target.value) || 30))
            displaydensity(ind)
            });
            function displaydensity(index) {
              if (!Number.isInteger(index) || !json_array[index] ||
                  json_array[index].length < 2) {
                Plotly.purge('secPlot');
                return;
              }
              var temp_x = json_array[index].map(coord => coord[0]);
              var temp_y = json_array[index].map(coord => coord[1]);
              var n = space;
              var shape = [n, 2];
              var strides = [1, n];
              var offset = 0;
              var order = 'column-major';

              var out = kde2d(temp_x, temp_y, {
                'n': n,
                'buffer': temp_x.concat(temp_y),
                'order': 'column-major',
                'offset': offset,
                'strides': strides
              });

              let twoDArray = [];
              for (let i = 0; i < out.z._buffer.length; i += n) {
                twoDArray.push(out.z._buffer.slice(i, i + n));
              }
              var arr = twoDArray;
              console.log(twoDArray)
              var data_u = [{
                z: twoDArray,
                type: 'surface'
              }];
              var layout_U = {
                title: 'kde2d plot for m/z: ' + data[index].mz,
                autosize: true,
                plot_bgcolor: 'rgba(0,0,0,0)', // Transparent background for the plot area
                paper_bgcolor: 'rgba(0,0,0,0)',
              };
              //loadingMessageElement.innerText = '';
              // Create dot plot
              Plotly.newPlot('secPlot', data_u, layout_U);
            }
            displaydensity(0)



            //table plot
            function displayTable(index) {
              var rowData = annotated_table[index];
              var tableData = [
                {
                  type: 'table',
                  header: {
                    values: ['mz', 'mz_name', 'annotated metabolites', 'entry of the metabolites'],
                    align: 'center',
                    line: {width: 1, color: 'black'},
                    fill: {color: 'grey'},
                    font: {family: 'Arial', size: 12, color: 'white'}
                  },
                  cells: {
                    values: rowData,
                    align: 'center',
                    line: {color: 'black', width: 1},
                    fill: {color: ['white', 'lightgrey']},
                    font: {family: 'Arial', size: 11, color: ['black']}
                  }
                }
              ];
              var layout_table = {
                title: 'Annotation table for m/z: ' + data[index].mz,
                autosize: true,
                plot_bgcolor: 'rgba(0,0,0,0)', // Transparent background for the plot area
                paper_bgcolor: 'rgba(0,0,0,0)',
              };

              Plotly.newPlot('tablePlot', tableData,layout_table);
            }
            displayTable(0);
            </script>
              </body>
              </html>"
  )

  if (!dir.exists(folder) && !dir.create(folder, recursive = TRUE))
    stop("Cannot create output directory.", call. = FALSE)
  output <- file.path(folder, "mzs_density_map.html")
  writeLines(javascript, output)
  invisible(output)
}






#' Interactively inspect mass windows in a spatial experiment
#'
#' The app sums every measured feature within the selected mass window without
#' changing the input assays or annotations. The mass-axis controls use rowData
#' masses, never mean intensities. No browser is launched by this function.
#' @param obj A SpatialExperiment or aligned Cardinal experiment.
#' @param assay Primary or alternative experiment.
#' @param slot Expression assay.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param image Optional image ID from imgData.
#' @param sampleId Optional sample ID. NULL displays separate sample panels.
#' @return A Shiny application object; use shiny::runApp to launch.
#' @rawNamespace import(shiny, except = runExample)
#' @export
#' @examples
#' utils::str(formals(interactiveSpatialPlot))
interactiveSpatialPlot <- function(obj, assay = "main", slot = "counts",
                                    image = NULL, sampleId = NULL) {
  obj <- if (is.null(sampleId)) .nativeSpatialObject(obj) else
    .nativeSingleSample(obj, sampleId)
  expression <- .nativeExpression(obj, assay, slot)
  masses <- .massValues(.featureMetadata(obj, assay), rownames(expression))
  if (!length(masses)) stop("No features to plot.", call. = FALSE)
  means <- as.numeric(Matrix::rowMeans(expression))
  if (any(!is.finite(means))) stop("Expression must be finite.", call. = FALSE)
  centre <- stats::median(masses)
  width <- max(diff(range(masses)) / 100, 0.001)
  ui <- shiny::fluidPage(
    shiny::titlePanel("Interactive Spatial Intensity Plot"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::numericInput("curve_center", "m/z value", value = centre,
          min = min(masses), max = max(masses)),
        shiny::radioButtons("bin_mode", "Window units",
          choices = c("Absolute (m/z)" = "mz", "Relative (ppm)" = "ppm"),
          selected = "mz"),
        shiny::numericInput("curve_width", "Half-width", value = width, min = 0),
        shiny::numericInput("spot_size", "Point size", value = 1.6, min = 0),
        shiny::numericInput("x_min", "Minimum m/z", value = min(masses)),
        shiny::numericInput("x_max", "Maximum m/z", value = max(masses))),
      shiny::mainPanel(
        plotly::plotlyOutput("plot"),
        shiny::verbatimTextOutput("selected_indices"),
        shiny::plotOutput("spatial_plot"))))
  server <- function(input, output, session) {
    window <- shiny::reactive({
      shiny::req(input$curve_center, !is.null(input$curve_width), input$bin_mode)
      .mzWindow(input$curve_center, input$curve_width, input$bin_mode)
    })
    selectedRows <- shiny::reactive({
      bounds <- window()
      which(masses >= bounds[1L] & masses <= bounds[2L])
    })
    selectedValues <- shiny::reactive({
      as.numeric(Matrix::colSums(expression[selectedRows(), , drop = FALSE]))
    })
    output$plot <- plotly::renderPlotly({
      bounds <- window()
      plot <- plotly::plot_ly(x = masses, y = means, type = "bar", name = "Mean intensity")
      plotly::layout(plot,
        xaxis = list(title = "m/z", range = c(input$x_min, input$x_max)),
        yaxis = list(title = "Mean intensity"),
        shapes = list(list(type = "rect", x0 = bounds[1L], x1 = bounds[2L],
          y0 = 0, y1 = 1, yref = "paper", fillcolor = "lightblue",
          opacity = 0.3, line = list(width = 0))))
    })
    output$selected_indices <- shiny::renderText({
      if (!length(selectedRows())) return("No features in this mass window.")
      paste("Selected features:", paste(rownames(expression)[selectedRows()], collapse = ", "))
    })
    output$spatial_plot <- shiny::renderPlot({
      shiny::validate(shiny::need(length(selectedRows()) > 0, "No selected features."))
      values <- matrix(selectedValues(), nrow = 1L,
        dimnames = list(paste0("m/z ", input$curve_center, " window"), colnames(obj)))
      .plotNativeExpression(obj, values, images = image, pointSize = input$spot_size)
    })
  }
  shiny::shinyApp(ui, server)
}

.mzWindow <- function(centre, width, units) {
  if (length(centre) != 1L || !is.finite(centre) || length(width) != 1L ||
      !is.finite(width) || width < 0 || !units %in% c("mz", "ppm"))
    stop("Supply a finite m/z centre and non-negative window width.", call. = FALSE)
  if (units == "ppm") {
    if (centre <= 0) stop("ppm windows require a positive mass.", call. = FALSE)
    width <- width * centre * 1e-6
  }
  c(centre - width, centre + width)
}
