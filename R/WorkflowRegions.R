.wf_region_summary <- function(E, labels, min_size) {
  labels <- as.character(labels)
  labelled <- !is.na(labels) & nzchar(trimws(labels))
  levels <- unique(labels[labelled])
  means <- detected <- matrix(NA_real_, nrow(E), length(levels),
    dimnames = list(rownames(E), levels))
  counts <- vapply(levels, function(g) sum(labelled & labels == g, na.rm = TRUE), integer(1))
  tables <- list()
  total <- as.numeric(Matrix::rowSums(E[, labelled, drop = FALSE]))
  positives <- as.numeric(Matrix::rowSums(E[, labelled, drop = FALSE] != 0))
  for (i in seq_along(levels)) {
    inside <- which(labelled & labels == levels[i])
    means[, i] <- as.numeric(Matrix::rowMeans(E[, inside, drop = FALSE]))
    detected[, i] <- as.numeric(Matrix::rowMeans(E[, inside, drop = FALSE] != 0))
    outside <- sum(labelled) - counts[i]
    if (counts[i] < min_size || outside < min_size) next
    reference <- (total - means[, i] * counts[i]) / outside
    reference_detected <- (positives - detected[, i] * counts[i]) / outside
    tables[[levels[i]]] <- data.frame(feature = rownames(E),
      effect = means[, i] - reference, mean_region = means[, i], mean_reference = reference,
      detected_region = detected[, i], detected_reference = reference_detected,
      n_region = counts[i], n_reference = outside, stringsAsFactors = FALSE)
  }
  list(field_levels = levels, means = means, detected = detected, tables = tables,
    groups = data.frame(region = levels, observations = unname(counts),
      eligible = counts >= min_size & sum(labelled) - counts >= min_size),
    unlabelled = sum(!labelled), method = "Descriptive region minus other labelled observations; no pixel-based P values. Non-zero fraction is measured on the workflow layer.")
}

.wf_preview_indices <- function(metadata, regions, maximum) {
  strata <- as.character(metadata$sample_id)
  if (!is.null(regions)) {
    label <- as.character(metadata[[regions]])
    label[is.na(label) | !nzchar(trimws(label))] <- "(unlabelled)"
    strata <- interaction(strata, label, drop = TRUE, lex.order = TRUE)
  }
  groups <- split(seq_len(nrow(metadata)), strata, drop = TRUE)
  # Allocate one point per stratum before adding further evenly spaced points.
  counts <- lengths(groups)
  quota <- pmin(counts, floor(maximum / max(1L, length(groups))))
  remaining <- min(sum(counts), maximum) - sum(quota)
  while (remaining > 0) {
    available <- head(which(quota < counts), remaining)
    quota[available] <- quota[available] + 1L
    remaining <- remaining - length(available)
  }
  sort(as.integer(unlist(lapply(seq_along(groups), function(i) {
    g <- groups[[i]]
    if (!quota[i]) return(integer())
    g[unique(as.integer(round(seq(1, length(g), length.out = quota[i]))))]
  }))))
}

.wf_modality_spatial <- function(object, assay, features, columns = seq_len(ncol(object))) {
  e <- .experimentForAssay(object, assay)
  SpatialExperiment::SpatialExperiment(
    assays = list(workflow = .assayData(object, assay, "workflow")[features, columns, drop = FALSE]),
    rowData = SummarizedExperiment::rowData(e)[features, , drop = FALSE],
    colData = SummarizedExperiment::colData(object)[columns, , drop = FALSE],
    spatialCoords = SpatialExperiment::spatialCoords(object)[columns, , drop = FALSE],
    metadata = S4Vectors::metadata(e))
}

.wf_native_analysis <- function(result, config) {
  defaults <- list(moran = FALSE, graph_pca = FALSE, max_points = 600L,
    max_features = 50L, neighbors = 4L, lambda = 0.5)
  if (!is.list(config) || (length(config) && (is.null(names(config)) ||
      any(!nzchar(names(config))) || anyDuplicated(names(config)))) ||
      length(setdiff(names(config), names(defaults))))
    stop("Unknown native analysis settings.", call. = FALSE)
  config <- utils::modifyList(defaults, config)
  for (n in c("moran", "graph_pca")) if (!is.logical(config[[n]]) || length(config[[n]]) != 1 || is.na(config[[n]]))
    stop(n, " must be TRUE or FALSE.", call. = FALSE)
  .wf_number(config$max_points, "native max_points", 4, TRUE)
  .wf_number(config$max_features, "native max_features", 2, TRUE)
  .wf_number(config$neighbors, "native neighbors", 1, TRUE)
  .wf_number(config$lambda, "native lambda", 0)
  xy <- SpatialExperiment::spatialCoords(result$object)
  answer <- list(config = config, modalities = list(), status = "completed")
  if (!config$moran && !config$graph_pca) {
    answer$status <- "not_requested"; answer$reason <- "Both native spatial modules are disabled."
    return(answer)
  }
  if (ncol(xy) < 2 || any(!is.finite(xy))) {
    answer$status <- "skipped"; answer$reason <- "No finite spatial coordinates supplied."
    return(answer)
  }
  sample_ids <- as.character(result$object$sample_id)
  for (name in names(result$analysis)) {
    features <- head(result$analysis[[name]]$pca$features, config$max_features)
    if (length(features) < 2) {
      answer$modalities[[name]] <- list(status = "skipped", reason = "Fewer than two variable features.")
      next
    }
    samples <- list()
    for (s in unique(sample_ids)) {
      cols <- which(sample_ids == s)
      cols <- cols[unique(as.integer(round(seq(1, length(cols), length.out = min(length(cols), config$max_points)))))]
      if (length(cols) < 4 || anyDuplicated(as.data.frame(xy[cols, , drop = FALSE]))) {
        samples[[s]] <- list(status = "skipped", reason = "Need four distinct spatial positions.")
        next
      }
      small <- .wf_modality_spatial(result$object, result$settings$assay_names[[name]], features, cols)
      info <- list(status = "completed", pixels = colnames(small), features = features,
        source_observations = sum(sample_ids == s), source_features = nrow(.experimentForAssay(result$object, result$settings$assay_names[[name]])))
      if (config$moran) {
        fit <- findSpatiallyVariableMetabolites(small, slot = "workflow", max_spots = NULL,
          nfeatures = length(features), verbose = FALSE)
        rd <- as.data.frame(SummarizedExperiment::rowData(fit))
        info$moran <- data.frame(feature = rownames(fit),
          I = rd$MoransI_observed, P.Value = rd$MoransI_p.value,
          FDR = rd$MoransI_p.adjust, sample_id = s, n_points = ncol(fit), family_size = nrow(fit))
        info$moran_provenance <- S4Vectors::metadata(fit)$moransi
      }
      if (config$graph_pca) {
        E <- t(scale(t(as.matrix(.assayData(small, layer = "workflow"))))); E[!is.finite(E)] <- 0
        SummarizedExperiment::assay(small, "scaled") <- E
        fit <- runSpatialGraphPCA(small, n_components = min(5L, nrow(E) - 1L),
          slot = "scaled", platform = "ST", n_neighbors = min(config$neighbors, ncol(E) - 1L),
          lambda = config$lambda, verbose = FALSE)
        info$graph_pca <- list(scores = SingleCellExperiment::reducedDim(fit, "SpatialPCA"),
          edges = SingleCellExperiment::colPair(fit, "SpatialKNN"),
          loadings = S4Vectors::metadata(fit)$reductionLoadings$SpatialPCA)
      }
      samples[[s]] <- info
    }
    answer$modalities[[name]] <- list(status = if (any(vapply(samples,
      function(s) identical(s$status, "completed"), logical(1)))) "completed" else "skipped", samples = samples)
  }
  if (!any(vapply(answer$modalities, function(m) identical(m$status, "completed"), logical(1)))) {
    answer$status <- "skipped"; answer$reason <- "No modality/sample had enough variable features and distinct spatial positions."
  }
  answer
}

#' Add region summaries and native spatial analyses to a workflow result
#'
#' Adds descriptive region-versus-rest effects, detection fractions, region
#' means and eligibility audits. Existing biological-replicate contrasts remain
#' unchanged. Optional native spatial modules call findSpatiallyVariableMetabolites
#' and runSpatialGraphPCA on an explicitly bounded, recorded feature/point screen.
#' Works with any modality names, species and observation metadata, including
#' results saved by earlier workflow versions.
#' @param result A spamtp_workflow result.
#' @param regions A colData column for regions; NULL uses the workflow group,
#'   or workflow_cluster if clustering was requested. No labels are invented.
#' @param min_region_size Minimum observations in both a region and its reference.
#' @param native Named settings: moran, graph_pca, max_points (600),
#'   max_features (50), neighbors (4), lambda (0.5). These additional bounded
#'   screens default FALSE; the workflow structure module runs full-observation
#'   Graph PCA separately and uses it in the comparative analysis.
#' @param markers Compute native scran marker scores through findAllDEMs.
#' @param spatial_blocks Number of coordinate-only blocks per sample for marker
#'   sensitivity; zero disables block omission. These are not biological replicates.
#' @return The result with region_analysis, native_analysis and updated provenance.
#' @details Region effects use all labelled retained observations. Blank/NA
#'   labels are excluded and counted in the audit. No P values are inferred from
#'   pixels or random pools. Moran's I statistics are exploratory and use BH
#'   within the selected feature family per modality and sample; this screening
#'   is not a genome-wide or biological-replicate test. Marker reference means
#'   equally weight other eligible regions, following scran pairwise summaries.
#' @export
#' @examples
#' utils::str(formals(analyzeSpaMTPRegions))
analyzeSpaMTPRegions <- function(result, regions = NULL, min_region_size = 3, native = list(),
    markers = TRUE, spatial_blocks = 6) {
  if (!inherits(result, "spamtp_workflow")) stop("result must be a spamtp_workflow.", call. = FALSE)
  .wf_number(min_region_size, "min_region_size", 1, TRUE)
  withr::local_seed(result$settings$seed)
  regions <- regions %||% result$settings$group
  if (is.null(regions) && "workflow_cluster" %in% colnames(SummarizedExperiment::colData(result$object)))
    regions <- "workflow_cluster"
  if (!is.null(regions) && (length(regions) != 1L || !regions %in% colnames(SummarizedExperiment::colData(result$object))))
    stop("regions must name one observation metadata field.", call. = FALSE)
  summaries <- list()
  if (!is.null(regions)) for (name in names(result$analysis)) {
    E <- .assayData(result$object, result$settings$assay_names[[name]], "workflow")
    summaries[[name]] <- .wf_region_summary(E, result$object[[regions]], min_region_size)
    if (isTRUE(markers)) {
      m <- findAllDEMs(result$object, ident = regions,
        assay = result$settings$assay_names[[name]], slot = "workflow", method = "markers",
        min_region_size = min_region_size, spatial_blocks = spatial_blocks, seed = result$settings$seed)
      summaries[[name]]$markers <- m
      if (identical(m$status, "completed")) {
        # Shared values feed ranking, spatial browsing and the native heatmap.
        summaries[[name]]$tables <- lapply(split(m$DEMs, m$DEMs$cluster), function(t) {
          t$feature <- t$gene; t$effect <- t$logFC
          t[order(-t$mean_auc, t$marker_rank, t$gene), , drop = FALSE]
        })
        summaries[[name]]$method <- m$method
      }
    }
  }
  result$region_analysis <- list(field = regions, min_region_size = min_region_size,
    modalities = summaries, status = if (is.null(regions)) "skipped" else "completed")
  result$native_analysis <- .wf_native_analysis(result, native)
  result$region_analysis$package_version <- result$native_analysis$package_version <-
    as.character(utils::packageVersion("SpaMTP"))
  result$region_analysis$created <- result$native_analysis$created <-
    format(Sys.time(), tz = "UTC", usetz = TRUE)
  result$settings$regions <- regions
  result$settings$markers <- list(enabled = markers, spatial_blocks = spatial_blocks,
    min_region_size = min_region_size)
  result$settings$native <- result$native_analysis$config
  result$stages <- result$stages[!result$stages$stage %in% c("region_summary", "native_spatial"), ]
  result$stages <- rbind(result$stages,
    data.frame(stage = "region_summary", status = result$region_analysis$status,
      detail = if (is.null(regions)) "No region metadata or requested clusters available." else paste("Region effects and detection on", regions, "; inferential tests retain their biological-replicate design.")),
    data.frame(stage = "native_spatial", status = result$native_analysis$status,
      detail = result$native_analysis$reason %||% "Native Moran's I and graph PCA settings, tested features and spatial subsets retained."))
  S4Vectors::metadata(result$object)$spamtp_workflow$stages <- result$stages
  S4Vectors::metadata(result$object)$spamtp_workflow$settings <- result$settings
  result
}
