.pathwayScores <- function(pathways, object, assay = "main", slot = "logcounts",
                           standardize = TRUE, pathway_index = NULL, duplicate_genes = "error") {
  expression <- .assayData(object, assay, slot)
  if (!is.list(pathways) || !length(pathways))
    stop("pathways must be a non-empty list of feature sets.", call. = FALSE)
  audit <- NULL
  if (!is.null(pathway_index)) {
    .pathway_validate_index(pathway_index)
    if (!is.null(pathway_index$gene_index)) .gene_check_experiment_species(object, assay, pathway_index$gene_index$provenance$organism)
    mapped <- .pathway_expression(expression, pathway_index, duplicate_genes)
    expression <- mapped$expression
    pathways <- lapply(pathways, function(ids) unique(c(ids[grepl("^RAMP_C_", ids)],
      .pathway_gene_ids(ids[!grepl("^RAMP_C_", ids)], pathway_index))))
    audit <- list(provenance = pathway_index$provenance, inputs = mapped$mapping,
      used_members = lapply(pathways, intersect, y = rownames(expression)))
  }
  scores <- lapply(pathways, function(features) {
    matched <- intersect(unique(features), rownames(expression))
    if (!length(matched))
      stop("A pathway has no features in the selected assay.", call. = FALSE)
    score <- as.numeric(Matrix::colSums(expression[matched, , drop = FALSE])) /
      sqrt(length(matched))
    if (any(!is.finite(score)))
      stop("Pathway scores require finite expression values.", call. = FALSE)
    score <- score - mean(score)
    if (isTRUE(standardize)) {
      deviation <- stats::sd(score)
      score <- if (is.finite(deviation) && deviation > 0) score / deviation else
        rep(0, length(score))
    }
    score
  })
  result <- do.call(cbind, scores)
  rownames(result) <- colnames(object)
  colnames(result) <- names(pathways) %||% paste0("pathway", seq_along(pathways))
  attr(result, "pathway_mapping") <- audit
  result
}

addGesecaScores <- function(pathways, object, assay = "main", slot = "logcounts",
                            prefix = "", scale = FALSE, pathway_index = NULL, duplicate_genes = "error") {
  object <- .nativeSpatialObject(object)
  scores <- .pathwayScores(pathways, object, assay, slot, standardize = scale,
    pathway_index = pathway_index, duplicate_genes = duplicate_genes)
  colnames(scores) <- paste0(prefix, colnames(scores))
  object <- .setCellMetadata(object, scores)
  S4Vectors::metadata(object)$pathway_score_mapping <- attr(scores, "pathway_mapping")
  object
}

.namedPathways <- function(pathways, database, database_version,
                           database_source, database_local_dir) {
  resources <- .spamtp_db_bundle(
    c("analytehaspathway", "pathway"), database = database,
    version = database_version, source = database_source,
    local_dir = database_local_dir)
  table <- merge(resources$analytehaspathway, resources$pathway,
                 by = "pathwayRampId")
  sets <- split(table$rampId, table$pathwayName)
  matched <- match(tolower(pathways), tolower(names(sets)))
  if (!length(pathways) || anyNA(matched))
    stop("Unknown pathway names: ", paste(pathways[is.na(matched)], collapse = ", "),
         call. = FALSE)
  stats::setNames(sets[matched], pathways)
}

#' Plot a pathway score on a Bioconductor dimensionality reduction
#'
#' Scores are the sum of matched feature values divided by the square root of
#' the matched feature count, then centred and scaled across pixels. Constant
#' scores are shown as zero. Plotting uses scater::plotReducedDim.
#' @inheritParams createPathwayAssay
#' @param pathway Character vector of feature IDs in one pathway.
#' @param object A SpatialExperiment with reducedDims.
#' @param title Optional plot title.
#' @param assay Primary (`"main"`) or alternative experiment name.
#' @param slot Expression assay name within the selected experiment
#'   (default = `"logcounts"`), read with `SummarizedExperiment::assay()`.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param reduction Name in reducedDimNames(object). NULL prefers UMAP, TSNE,
#'   then PCA, otherwise the first available reduction.
#' @param colors Gradient colours.
#' @param guide Colour guide type.
#' @param ... Additional arguments to scater::plotReducedDim.
#' @return A ggplot displaying the pathway score.
#' @export
#' @examples
#' x <- SpatialExperiment::SpatialExperiment(
#'   assays = list(logcounts = rbind(a = 1:4, b = 4:1)),
#'   spatialCoords = cbind(x = 1:4, y = 0),
#'   colData = S4Vectors::DataFrame(row.names = paste0("p", 1:4)))
#' SingleCellExperiment::reducedDim(x, "PCA") <- cbind(1:4, c(1, 3, 2, 4))
#' plotSinglePathway("a", x)
plotSinglePathway <- function(pathway, object, title = NULL, assay = "main",
                              slot = "logcounts", reduction = NULL,
                              colors = c("darkblue", "lightgrey", "darkred"),
                              guide = "colourbar", pathway_index = NULL, duplicate_genes = "error", ...) {
  object <- .nativeSpatialObject(object)
  scores <- .pathwayScores(list(pathway = pathway), object, assay, slot,
    pathway_index = pathway_index, duplicate_genes = duplicate_genes)
  plot <- .plot_pathway_reduced_score(object, scores[, 1L], title, reduction, colors, guide, ...)
  attr(plot, "pathway_mapping") <- attr(scores, "pathway_mapping")
  plot
}

.plot_pathway_reduced_score <- function(object, score, title, reduction, colors, guide, ...) {
  column <- tail(make.unique(c(colnames(SummarizedExperiment::colData(object)),
                               rownames(object), "pathway_score")), 1L)
  SummarizedExperiment::colData(object)[[column]] <- score
  available <- SingleCellExperiment::reducedDimNames(object)
  if (is.null(reduction) && length(available)) {
    preferred <- match(c("umap", "tsne", "pca"), tolower(available))
    preferred <- preferred[!is.na(preferred)]
    reduction <- available[if (length(preferred)) preferred[[1L]] else 1L]
  }
  if (length(reduction) != 1L || !reduction %in% available)
    stop("Choose a reduction from reducedDimNames(object).", call. = FALSE)
  suppressMessages(scater::plotReducedDim(
    object, dimred = reduction, colour_by = column, ...) +
      ggplot2::scale_colour_gradientn(colors = colors, limits = c(-3, 3),
        oob = scales::squish, guide = guide, name = "z-score") +
      ggplot2::labs(title = title))
}

#' Plot named pathways on a dimensionality reduction
#' @param pathways Character vector of pathwayRampId values or unambiguous database pathway names.
#' @inheritParams plotSinglePathway
#' @inheritParams createPathwayAssay
#' @param database Optional named list of database tables.
#' @param database_version Database version.
#' @param database_source Database source, auto or spamtpdb.
#' @param database_local_dir Optional local database directory.
#' @return A named list of ggplots with pathway_coverage and pathway_index attributes; each plot also carries its exact pathway_scores.
#' @export
plotPathways <- function(pathways, object, title = NULL, assay = "main",
    slot = "logcounts", reduction = NULL, colors = c("darkblue", "lightgrey", "darkred"),
    guide = "colourbar", database = NULL, database_version = "latest",
    database_source = c("auto", "spamtpdb"), database_local_dir = NULL,
    pathway_index = NULL, gene_mapping = c("auto", "hgnc", "ramp"), gene_reference = NULL,
    gene_index = NULL, gene_reference_version = "latest", gene_reference_local_dir = NULL,
    organism = "Homo sapiens", duplicate_genes = c("error", "mean", "sum"), ...) {
  object <- .nativeSpatialObject(object)
  context <- .pathway_score_context(object, assay, slot, database, pathway_index,
    match.arg(gene_mapping), gene_reference, gene_index, gene_reference_version,
    gene_reference_local_dir, organism, database_version, match.arg(database_source),
    database_local_dir, match.arg(duplicate_genes))
  ids <- .pathway_select(context$index, pathways)
  selected <- context$coverage[match(ids, context$coverage$pathwayRampId), , drop = FALSE]
  if (any(selected$used_size == 0)) stop("A requested pathway has no measured members.", call. = FALSE)
  scores <- .pathway_score_matrix(context$expression, .pathway_sets_from_coverage(selected))
  titles <- title %||% selected$pathwayName
  if (length(titles) != length(ids)) stop("Provide one title per pathway.", call. = FALSE)
  plots <- lapply(seq_along(ids), function(i) {
    p <- .plot_pathway_reduced_score(object, scores[i, ], titles[i], reduction, colors, guide, ...)
    attr(p, "pathway_coverage") <- selected[i, , drop = FALSE]
    attr(p, "pathway_scores") <- scores[i, ]
    p
  })
  names(plots) <- pathways
  attr(plots, "pathway_coverage") <- selected
  attr(plots, "pathway_index") <- context$index$provenance
  plots
}
#' Plot a pathway score in tissue coordinates
#' @inheritParams plotSinglePathway
#' @param images Image IDs in imgData to overlay; NULL draws coordinates only.
#' @param image.alpha Image opacity.
#' @param crop Whether to restrict plotting to the pixel extent.
#' @param min.cutoff,max.cutoff Numeric or quantile colour cutoffs.
#' @param ncol Number of columns in the combined plot.
#' @param pt.size.factor Point size.
#' @param alpha Point opacity (the first value is used).
#' @param shape Point shape.
#' @param stroke Point border width.
#' @param interactive Convert the result to an interactive plotly plot.
#' @param image.labels Optional labels for the selected images.
#' @return A ggplot or combined plot; a plotly object if interactive is TRUE.
#' @export
#' @examples
#' x <- SpatialExperiment::SpatialExperiment(
#'   assays = list(logcounts = rbind(a = 1:4, b = 4:1)),
#'   spatialCoords = cbind(x = 1:4, y = 0),
#'   colData = S4Vectors::DataFrame(row.names = paste0("p", 1:4)))
#' plotSinglePathwaySpatially("a", x)
plotSinglePathwaySpatially <- function(
    pathway, object, images = NULL, title = NULL, image.alpha = 1,
    assay = "main", slot = "logcounts",
    colors = c("darkblue", "lightgrey", "darkred"), guide = "colourbar",
    crop = TRUE, min.cutoff = NA, max.cutoff = NA, ncol = NULL,
    pt.size.factor = 1.6, alpha = 1, shape = 16, stroke = 0,
    interactive = FALSE, image.labels = NULL, pathway_index = NULL, duplicate_genes = "error"
) {
  object <- .nativeSpatialObject(object)
  values <- .pathwayScores(list(pathway = pathway), object, assay, slot,
    pathway_index = pathway_index, duplicate_genes = duplicate_genes)
  scores <- values[, 1L]
  plot <- .plot_pathway_spatial_score(scores, object, images, title, image.alpha,
    colors, guide, crop, min.cutoff, max.cutoff, ncol, pt.size.factor,
    alpha, shape, stroke, interactive, image.labels)
  attr(plot, "pathway_mapping") <- attr(values, "pathway_mapping")
  plot
}

.plot_pathway_spatial_score <- function(scores, object, images = NULL,
    title = NULL, image.alpha = 1,
    colors = c("darkblue", "lightgrey", "darkred"), guide = "colourbar",
    crop = TRUE, min.cutoff = NA, max.cutoff = NA, ncol = NULL,
    pt.size.factor = 1.6, alpha = 1, shape = 16, stroke = 0,
    interactive = FALSE, image.labels = NULL) {
  .plotNativeSpatialScore(
    object, scores, images = images, title = title, colors = colors,
    guide = guide, imageAlpha = image.alpha, pointSize = pt.size.factor,
    alpha = alpha[[1L]], shape = shape, stroke = stroke, crop = crop,
    minCutoff = min.cutoff, maxCutoff = max.cutoff, ncol = ncol,
    interactive = interactive, imageLabels = image.labels)
}

#' Plot named pathways in tissue coordinates
#' @inheritParams plotSinglePathwaySpatially
#' @inheritParams plotPathways
#' @param ... Additional arguments to plotSinglePathwaySpatially.
#' @return A named list of spatial pathway plots, carrying the same coverage and score attributes as plotPathways().
#' @export
plotPathwaysSpatially <- function(pathways, object, images = NULL, title = NULL,
    image.alpha = 1, assay = "main", slot = "logcounts",
    colors = c("darkblue", "lightgrey", "darkred"), guide = "colourbar",
    database = NULL, database_version = "latest", database_source = c("auto", "spamtpdb"),
    database_local_dir = NULL, pathway_index = NULL, gene_mapping = c("auto", "hgnc", "ramp"),
    gene_reference = NULL, gene_index = NULL, gene_reference_version = "latest",
    gene_reference_local_dir = NULL, organism = "Homo sapiens",
    duplicate_genes = c("error", "mean", "sum"), ...) {
  object <- .nativeSpatialObject(object)
  context <- .pathway_score_context(object, assay, slot, database, pathway_index,
    match.arg(gene_mapping), gene_reference, gene_index, gene_reference_version,
    gene_reference_local_dir, organism, database_version, match.arg(database_source),
    database_local_dir, match.arg(duplicate_genes))
  ids <- .pathway_select(context$index, pathways)
  selected <- context$coverage[match(ids, context$coverage$pathwayRampId), , drop = FALSE]
  if (any(selected$used_size == 0)) stop("A requested pathway has no measured members.", call. = FALSE)
  scores <- .pathway_score_matrix(context$expression, .pathway_sets_from_coverage(selected))
  titles <- title %||% selected$pathwayName
  if (length(titles) != length(ids)) stop("Provide one title per pathway.", call. = FALSE)
  plots <- lapply(seq_along(ids), function(i) {
    p <- .plot_pathway_spatial_score(scores[i, ], object, images, titles[i], image.alpha,
      colors, guide, ...)
    attr(p, "pathway_coverage") <- selected[i, , drop = FALSE]
    attr(p, "pathway_scores") <- scores[i, ]
    p
  })
  names(plots) <- pathways
  attr(plots, "pathway_coverage") <- selected
  attr(plots, "pathway_index") <- context$index$provenance
  plots
}
