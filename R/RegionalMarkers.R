.replicateRegionContrasts <- function(E, metadata, group, replicate, contrasts) {
  if (is.null(group)) return(list())
  if (!group %in% names(metadata) || anyNA(metadata[[group]]))
    stop("group must name complete observation metadata.", call. = FALSE)
  labels <- as.character(metadata[[group]])
  means <- matrix(vapply(unique(labels), function(g)
    as.numeric(Matrix::rowMeans(E[, labels == g, drop = FALSE])),
    numeric(nrow(E))), nrow = nrow(E), ncol = length(unique(labels)))
  dimnames(means) <- list(rownames(E), unique(labels))
  out <- list(group_means = means, tests = list())
  if (!length(contrasts)) return(out)
  if (is.null(replicate) || !replicate %in% names(metadata) || anyNA(metadata[[replicate]]))
    stop("Inferential contrasts require a complete biological replicate column.", call. = FALSE)
  ids <- as.character(metadata[[replicate]])
  .wf_named(contrasts, "contrasts")
  for (name in names(contrasts)) {
    c <- contrasts[[name]]
    if (length(setdiff(names(c), c("numerator", "denominator", "paired"))) ||
        length(c$numerator) != 1L || length(c$denominator) != 1L ||
        c$numerator == c$denominator || !all(c(c$numerator, c$denominator) %in% labels))
      stop("Each contrast needs distinct observed numerator and denominator groups.", call. = FALSE)
    keep <- labels %in% c(c$numerator, c$denominator)
    units <- unique(data.frame(replicate = ids[keep], group = labels[keep]))
    paired <- isTRUE(c$paired)
    if (paired) {
      complete <- intersect(units$replicate[units$group == c$numerator],
        units$replicate[units$group == c$denominator])
      units <- units[units$replicate %in% complete, , drop = FALSE]
    } else if (anyDuplicated(units$replicate)) {
      stop("Shared biological replicates require paired=TRUE: ", name, call. = FALSE)
    }
    counts <- table(factor(units$group, levels = c(c$numerator, c$denominator)))
    if (any(counts < 3L)) {
      out$tests[[name]] <- list(status = "skipped",
        reason = "At least three biological replicates per group (three complete pairs if paired) are required.",
        units = units, contrast = c)
      next
    }
    Y <- matrix(vapply(seq_len(nrow(units)), function(i) as.numeric(Matrix::rowMeans(
      E[, ids == units$replicate[i] & labels == units$group[i], drop = FALSE])),
      numeric(nrow(E))), nrow = nrow(E), ncol = nrow(units))
    design <- cbind(denominator = as.numeric(units$group == c$denominator),
      numerator = as.numeric(units$group == c$numerator))
    if (paired) design <- cbind(design, stats::model.matrix(~ factor(units$replicate))[, -1, drop = FALSE])
    if (qr(design)$rank < ncol(design) || nrow(design) <= ncol(design))
      stop("Contrast has a rank-deficient design or no residual degrees of freedom.", call. = FALSE)
    fit <- limma::lmFit(Y, design)
    fit <- limma::contrasts.fit(fit, c(-1, 1, rep(0, ncol(design) - 2)))
    fit <- limma::eBayes(fit, trend = TRUE)
    tab <- limma::topTable(fit, number = Inf, sort.by = "none", confint = TRUE)
    tab$feature <- rownames(E)
    rownames(Y) <- rownames(E)
    colnames(Y) <- paste(units$replicate, units$group, sep = "::")
    out$tests[[name]] <- list(status = "completed", table = tab, units = units,
      expression = Y, design = design, coefficient = c(-1, 1, rep(0, ncol(design) - 2)),
      residual_df = fit$df.residual,
      contrast = c, method = "limma on equally weighted biological-replicate means of workflow values; BH within modality and contrast")
  }
  out
}

.nativeSpatialBlocks <- function(data, blocks, seed) {
  answer <- rep(NA_character_, ncol(data))
  if (!blocks) return(answer)
  if (!methods::is(data, "SpatialExperiment")) return(answer)
  xy <- SpatialExperiment::spatialCoords(data)
  if (ncol(xy) < 2L || any(!is.finite(xy))) return(answer)
  samples <- as.character(data$sample_id)
  withr::local_seed(seed)
  for (s in unique(samples)) {
    i <- which(samples == s)
    k <- min(blocks, nrow(unique(xy[i, , drop = FALSE])), floor(length(i) / 3))
    if (k < 2L) next
    z <- stats::kmeans(xy[i, , drop = FALSE], centers = k, nstart = 10, iter.max = 100)$cluster
    answer[i] <- paste(s, z, sep = "::")
  }
  answer
}

.nativeRegionMarkers <- function(data, ident, assay, slot, min_region_size,
    spatial_blocks, seed, annotation.column) {
  .validateFeatureCount(min_region_size, "min_region_size")
  if (!is.numeric(spatial_blocks) || length(spatial_blocks) != 1L ||
      is.na(spatial_blocks) || spatial_blocks < 0 || spatial_blocks != floor(spatial_blocks))
    stop("spatial_blocks must be a non-negative integer.", call. = FALSE)
  E <- .assayData(data, assay, slot)
  if (is.null(rownames(E)) || anyDuplicated(rownames(E)) || any(!is.finite(E)))
    stop("Marker analysis requires unique feature IDs and finite values.", call. = FALSE)
  labels <- as.character(SummarizedExperiment::colData(data)[[ident]])
  valid <- !is.na(labels) & nzchar(trimws(labels))
  counts <- table(labels[valid])
  groups <- names(counts)[counts >= min_region_size]
  valid <- valid & labels %in% groups
  eligibility <- data.frame(region = names(counts), observations = as.integer(counts),
    eligible = names(counts) %in% groups)
  if (length(groups) < 2L) return(list(status = "skipped", groups = eligibility,
    reason = "Need at least two regions with min_region_size observations."))
  upstream_deprecations <- character()
  score <- function(E, groups) withCallingHandlers(scran::scoreMarkers(E, groups = groups),
    warning = function(w) {
      text <- conditionMessage(w)
      if (grepl("^'(computeMinRank|summarizeAssayByGroup)' is deprecated", text)) {
        upstream_deprecations <<- unique(c(upstream_deprecations, text))
        invokeRestart("muffleWarning")
      }
    })
  fit <- score(E[, valid, drop = FALSE], groups = labels[valid])
  blocks <- .nativeSpatialBlocks(data, spatial_blocks, seed)
  valid_blocks <- unique(stats::na.omit(blocks[valid]))
  omissions <- list()
  # The same coordinate-only partition is used for every modality and region.
  # Refit marker scores, conditional on the original region labels; no P values.
  for (b in valid_blocks) {
    keep <- valid & (is.na(blocks) | blocks != b)
    n <- table(factor(labels[keep], levels = groups))
    if (any(n < min_region_size)) next
    omissions[[b]] <- score(E[, keep, drop = FALSE], groups = labels[keep])
  }
  tables <- lapply(groups, function(g) {
    f <- as.data.frame(fit[[g]])
    self <- which(valid & labels == g)
    other <- setdiff(groups, g)
    other_detected <- rowMeans(matrix(vapply(other, function(o)
      as.numeric(Matrix::rowMeans(E[, valid & labels == o, drop = FALSE] != 0)),
      numeric(nrow(E))), nrow = nrow(E), ncol = length(other)))
    tab <- data.frame(gene = rownames(E), cluster = g,
      logFC = f$self.average - f$other.average,
      mean_region = f$self.average, mean_reference = f$other.average,
      detected_region = as.numeric(Matrix::rowMeans(E[, self, drop = FALSE] != 0)),
      detected_reference = other_detected,
      mean_cohen = f$mean.logFC.cohen, min_cohen = f$min.logFC.cohen,
      mean_auc = f$mean.AUC, min_auc = f$min.AUC,
      n_region = length(self), n_reference = sum(valid) - length(self))
    tab$marker_rank <- rank(-tab$mean_auc, ties.method = "min", na.last = "keep")
    tab$stability_fits <- length(omissions)
    tab$direction_stability <- tab$cohen_min_block <- tab$cohen_max_block <- NA_real_
    tab$rank_best_block <- tab$rank_worst_block <- NA_real_
    if (length(omissions)) {
      d <- matrix(vapply(omissions, function(z) as.numeric(z[[g]]$mean.logFC.cohen),
        numeric(nrow(E))), nrow = nrow(E), ncol = length(omissions))
      auc <- matrix(vapply(omissions, function(z) as.numeric(z[[g]]$mean.AUC),
        numeric(nrow(E))), nrow = nrow(E), ncol = length(omissions))
      ranks <- matrix(apply(auc, 2, function(z)
        rank(-z, ties.method = "min", na.last = "keep")),
        nrow = nrow(E), ncol = length(omissions))
      tab$direction_stability <- rowMeans(sign(d) == sign(tab$mean_cohen), na.rm = TRUE)
      tab$cohen_min_block <- apply(d, 1, min, na.rm = TRUE)
      tab$cohen_max_block <- apply(d, 1, max, na.rm = TRUE)
      tab$rank_best_block <- apply(ranks, 1, min, na.rm = TRUE)
      tab$rank_worst_block <- apply(ranks, 1, max, na.rm = TRUE)
    }
    tab
  })
  dems <- do.call(rbind, tables)
  if (!is.null(annotation.column)) {
    rd <- .featureMetadata(data, assay)
    if (!annotation.column %in% names(rd)) stop("Unknown annotation.column.", call. = FALSE)
    dems[[annotation.column]] <- rd[[annotation.column]][match(dems$gene, rownames(E))]
  }
  means <- matrix(vapply(groups, function(g)
    as.numeric(Matrix::rowMeans(E[, valid & labels == g, drop = FALSE])),
    numeric(nrow(E))), nrow = nrow(E), ncol = length(groups))
  dimnames(means) <- list(rownames(E), groups)
  structure(list(status = "completed", DEMs = dems, scores = fit,
    expression = means, samples = data.frame(ident = groups, row.names = groups),
    groups = eligibility, blocks = data.frame(observation = colnames(E), block = blocks),
    method = "scran::scoreMarkers: pairwise AUC and Cohen effect sizes, equally summarized across other eligible regions; no pixel P values",
    provenance = list(function_name = "findAllDEMs", method = "markers", ident = ident,
      assay = assay, layer = slot, features = rownames(E), observations = colnames(E)[valid],
      min_region_size = min_region_size, spatial_blocks = spatial_blocks,
      block_omissions = names(omissions), block_omissions_ineligible = setdiff(valid_blocks, names(omissions)),
      stability = "Leave one coordinate-only k-means block out; region labels held fixed. Sensitivity, not independent replication or confidence intervals.",
      seed = seed, upstream_deprecations = upstream_deprecations,
      scran_version = as.character(utils::packageVersion("scran")))),
    class = c("spamtp_markers", "list"))
}
