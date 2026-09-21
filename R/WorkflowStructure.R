.wf_structure_config <- function(config, clusters) {
  defaults <- list(spatial = TRUE, primary = "pca", neighbors = 6L, lambda = 0.5,
    platform = "Visium", reference = NULL, k_grid = clusters,
    stability_blocks = 6L, umap = TRUE)
  if (!is.list(config) || length(setdiff(names(config), names(defaults))))
    stop("Unknown structure settings.", call. = FALSE)
  config <- utils::modifyList(defaults, config)
  config$primary <- match.arg(config$primary, c("pca", "spatial"))
  config$platform <- match.arg(config$platform, c("Visium", "ST"))
  for (n in c("spatial", "umap")) if (!is.logical(config[[n]]) ||
    length(config[[n]]) != 1L || is.na(config[[n]])) stop("Invalid ", n, call. = FALSE)
  .wf_number(config$neighbors, "structure neighbors", 1, TRUE)
  .wf_number(config$lambda, "structure lambda", 0)
  .wf_number(config$stability_blocks, "stability_blocks", 0, TRUE)
  if (length(config$k_grid) && (!is.numeric(config$k_grid) || anyNA(config$k_grid) ||
    any(config$k_grid < 2 | config$k_grid != floor(config$k_grid)))) stop("Invalid k_grid.", call. = FALSE)
  config$k_grid <- sort(unique(c(config$k_grid, clusters)))
  config
}

.wf_native_cluster <- function(z, k, seed) {
  x <- SingleCellExperiment::SingleCellExperiment(
    colData = S4Vectors::DataFrame(row.names = rownames(z)),
    reducedDims = S4Vectors::SimpleList(embedding = z))
  x <- getKmeanClusters(x, reduction = "embedding", cluster.name = "cluster",
    clusters = k, nstart = 20, iter.max = 100, seed = seed)
  as.integer(x$cluster)
}

.wf_structure <- function(result, config) {
  config <- .wf_structure_config(config, result$settings$clusters)
  object <- result$object
  xy <- SpatialExperiment::spatialCoords(object)
  spatial <- ncol(xy) >= 2L && all(is.finite(xy))
  if (config$primary == "spatial" && (!spatial || !config$spatial))
    stop("A spatial primary representation requires coordinates and spatial=TRUE.", call. = FALSE)
  if (!is.null(config$reference) && (length(config$reference) != 1L ||
      !config$reference %in% colnames(SummarizedExperiment::colData(object))))
    stop("structure reference must name external annotation metadata.", call. = FALSE)
  reference <- if (is.null(config$reference)) NULL else as.character(object[[config$reference]])
  if (identical(config$reference, "workflow_cluster"))
    stop("The evaluated partition cannot also be its reference.", call. = FALSE)
  representations <- list(); graph_names <- list(); graph <- NULL
  for (name in names(result$analysis)) {
    p <- result$analysis[[name]]$pca
    if (!identical(p$status, "completed")) next
    representations[[paste(name, "PCA", sep = " / ")]] <- list(embedding = p$scores,
      modality = name, method = "PCA", features = p$features,
      reduction = paste0("workflow_", name))
    if (config$spatial && spatial) {
      small <- .wf_modality_spatial(object, result$settings$assay_names[[name]], p$features)
      small <- scaleSMData(small, slot = "workflow")
      small <- runSpatialGraphPCA(small, n_components = ncol(p$scores),
        platform = config$platform, n_neighbors = config$neighbors, lambda = config$lambda,
        verbose = FALSE)
      scores <- SingleCellExperiment::reducedDim(small, "SpatialPCA")
      reduction <- paste0("workflow_graph_", name)
      SingleCellExperiment::reducedDim(object, reduction) <- scores
      assay <- result$settings$assay_names[[name]]
      if (assay != "main") {
        e <- SingleCellExperiment::altExp(object, assay)
        SingleCellExperiment::reducedDim(e, reduction) <- scores
        SingleCellExperiment::altExp(object, assay) <- e
      }
      graph_names[[name]] <- reduction
      result$analysis[[name]]$spatial_pca <- list(scores = scores,
        features = p$features, loadings = S4Vectors::metadata(small)$reductionLoadings$SpatialPCA,
        provenance = S4Vectors::metadata(small)$reduction_inputs$SpatialPCA)
      representations[[paste(name, "Graph PCA", sep = " / ")]] <- list(
        embedding = scores, modality = name, method = "Graph PCA", features = p$features, reduction = reduction)
      graph <- S4Vectors::metadata(small)$spatialGraphs$SpatialKNN
    }
  }
  if (!is.null(result$joint)) representations[["Joint / PCA"]] <- list(
    embedding = result$joint$embedding, modality = "Joint", method = "PCA", reduction = "integrated")
  if (length(graph_names) >= 2L && length(graph_names) == length(result$analysis)) {
    previous <- SingleCellExperiment::reducedDim(object, "integrated")
    object <- multiOmicIntegration(object, modalities = unname(result$settings$assay_names),
      reduction.list = unname(graph_names), dims.list = lapply(graph_names, function(n)
        seq_len(ncol(SingleCellExperiment::reducedDim(object, n)))), return.intermediate = TRUE)
    z <- SingleCellExperiment::reducedDim(object, "integrated")
    SingleCellExperiment::reducedDim(object, "workflow_joint_spatial") <- z
    SingleCellExperiment::reducedDim(object, "integrated") <- previous
    representations[["Joint / Graph PCA"]] <- list(embedding = z, modality = "Joint",
      method = "Graph PCA", reduction = "workflow_joint_spatial")
  }
  blocks <- .nativeSpatialBlocks(object, config$stability_blocks, result$settings$seed)
  omissions <- unique(stats::na.omit(blocks))
  metrics <- list(); assignments <- list(); stability <- list()
  withr::local_seed(result$settings$seed)
  for (name in names(representations)) {
    z <- as.matrix(representations[[name]]$embedding)
    # UMAP is a display of the entire representation; clustering never uses it.
    display <- if (ncol(z) == 1L) cbind(z[, 1], 0) else if (config$umap && nrow(z) >= 30L) scater::calculateUMAP(z,
      transposed = TRUE, n_neighbors = min(15L, nrow(z) - 1L), n_threads = 1L,
      n_sgd_threads = 1L, verbose = FALSE) else stats::prcomp(z, rank. = 2)$x[, 1:2, drop = FALSE]
    rownames(display) <- rownames(z)
    representations[[name]]$display <- display
    representations[[name]]$display_method <- if (config$umap && nrow(z) >= 30L && ncol(z) > 1L) "UMAP" else "PCA display"
    assignments[[name]] <- list()
    for (k in config$k_grid) {
      if (k >= nrow(unique(z))) next
      label <- .wf_native_cluster(z, k, result$settings$seed)
      assignments[[name]][[as.character(k)]] <- label
      silhouette <- bluster::approxSilhouette(z, factor(label))$width
      ari <- NA_real_
      if (!is.null(reference)) {
        valid <- !is.na(reference) & nzchar(trimws(reference))
        if (sum(valid) > 2L && length(unique(reference[valid])) > 1L)
          ari <- mclust::adjustedRandIndex(label[valid], reference[valid])
      }
      conditional <- numeric()
      for (b in omissions) {
        keep <- which(is.na(blocks) | blocks != b)
        if (k >= nrow(unique(z[keep, , drop = FALSE]))) next
        predicted <- .wf_native_cluster(z[keep, , drop = FALSE], k, result$settings$seed)
        conditional[b] <- mclust::adjustedRandIndex(label[keep], predicted)
      }
      stability[[paste(name, k, sep = "::")]] <- conditional
      spatial_agreement <- NA_real_
      if (!is.null(graph)) {
        edges <- summary(methods::as(graph, "generalMatrix"))
        if (nrow(edges)) spatial_agreement <- stats::weighted.mean(label[edges$i] == label[edges$j], edges$x)
      }
      metrics[[length(metrics) + 1L]] <- data.frame(representation = name, k = k,
        approx_silhouette = mean(silhouette), reference_ARI = ari,
        conditional_stability = if (length(conditional)) stats::median(conditional) else NA_real_,
        stability_omissions = length(conditional), spatial_neighbor_agreement = spatial_agreement,
        smallest_cluster = min(table(label)))
    }
  }
  primary <- paste(if (!is.null(result$joint)) "Joint" else names(result$analysis)[1],
    if (config$primary == "spatial") "Graph PCA" else "PCA", sep = " / ")
  if (!is.null(result$settings$clusters) && primary %in% names(assignments)) {
    label <- assignments[[primary]][[as.character(result$settings$clusters)]]
    object$workflow_cluster <- factor(label)
    result$clustering <- list(labels = label, centers = result$settings$clusters,
      input = primary, method = "getKmeanClusters on the full primary representation; k specified before evaluation")
    if (!is.null(result$joint)) result$joint$clusters <- label
  }
  if (!is.null(result$joint) && config$primary == "spatial") {
    result$joint$embedding <- representations[[primary]]$embedding
    result$joint$display <- representations[[primary]]$display
    result$joint$method <- "multiOmicIntegration of full-observation, feature-matched graph PCA embeddings; equal modality weights"
    SingleCellExperiment::reducedDim(object, "integrated") <- result$joint$embedding
  }
  result$object <- object
  result$structure <- list(config = config, representations = representations,
    metrics = if (length(metrics)) do.call(rbind, metrics) else data.frame(),
    assignments = assignments, stability = stability,
    blocks = data.frame(observation = colnames(object), block = blocks),
    primary = primary, graph = graph,
    interpretation = paste("Feature selection and scaling are identical for PCA and Graph PCA within each modality.",
      "Graph PCA uses all retained observations with separate spatial graphs per sample.",
      "UMAP is for display only. K and the primary method are not selected using reference ARI.",
      "Approximate silhouettes use centroid distances. Spatial agreement rewards smoothness by construction.",
      "Block omission refits clustering conditional on full-data embeddings; it does not validate preprocessing or generalization to new specimens."))
  result$settings$structure <- config
  result$stages <- rbind(result$stages, data.frame(stage = "structure_evaluation", status = "completed",
    detail = paste(length(representations), "full-observation representations; primary:", primary,
      "; reference is used only for evaluation; K is not selected from reference ARI.")))
  if (!is.null(result$joint)) result$stages$detail[result$stages$stage == "integration"] <- result$joint$method
  result
}
