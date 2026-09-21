.wf_explorer_payload <- function(result, points, features, rows) {
  metadata <- as.data.frame(SummarizedExperiment::colData(result$object))
  field <- result$region_analysis$field
  ids <- .wf_preview_indices(metadata, field, points)
  xy <- SpatialExperiment::spatialCoords(result$object)
  coords <- if (ncol(xy) >= 2L && all(is.finite(xy))) unname(xy[ids, 1:2, drop = FALSE]) else NULL
  vector <- function(x) unname(as.list(x))
  payload <- list(schema = 1L, observations = vector(colnames(result$object)[ids]),
    total_observations = ncol(result$object), coordinates = coords,
    samples = vector(as.character(metadata$sample_id[ids])),
    regions = if (is.null(field)) vector(rep("unlabelled", length(ids))) else vector(as.character(metadata[[field]][ids])),
    conditions = if (is.null(result$settings$group)) NULL else vector(as.character(metadata[[result$settings$group]][ids])),
    region_field = field, units = result$settings$coordinate_units, modalities = list())
  for (name in names(result$analysis)) {
    E <- .assayData(result$object, result$settings$assay_names[[name]], "workflow")
    rd <- .featureMetadata(result$object, result$settings$assay_names[[name]])
    labels <- rownames(E)
    label_column <- intersect(c("symbol", "gene_name", "all_IsomerNames"), names(rd))
    if (length(label_column)) {
      candidates <- as.character(rd[[label_column[1]]])
      valid <- !is.na(candidates) & nzchar(candidates) & !candidates %in% c("No Annotation", "NA")
      labels[valid] <- paste0(labels[valid], " | ", substr(candidates[valid], 1, 120))
    }
    region <- result$region_analysis$modalities[[name]]
    modes <- list(); priority <- list()
    add_mode <- function(table, title, kind, target, reference, method) {
      f <- match(table$feature, rownames(E))
      if (anyNA(f)) stop("Explorer features do not match the selected modality.", call. = FALSE)
      value <- function(n) if (n %in% names(table)) table[[n]] else rep(NA_real_, nrow(table))
      tab <- cbind(feature = f - 1L, effect = value("effect"),
        mean_region = value("mean_region"), mean_reference = value("mean_reference"),
        detected_region = value("detected_region"), detected_reference = value("detected_reference"),
        FDR = value("FDR"), P.Value = value("P.Value"), mean_auc = value("mean_auc"),
        min_auc = value("min_auc"), mean_cohen = value("mean_cohen"), min_cohen = value("min_cohen"),
        direction_stability = value("direction_stability"), rank_worst_block = value("rank_worst_block"),
        cohen_min_block = value("cohen_min_block"), cohen_max_block = value("cohen_max_block"))
      order <- if (kind == "test") order(tab[, "FDR"], -abs(tab[, "effect"]), na.last = TRUE) else
        if (any(is.finite(tab[, "mean_auc"]))) order(-tab[, "mean_auc"], -tab[, "effect"], na.last = TRUE) else order(-abs(tab[, "effect"]), na.last = TRUE)
      selected <- head(order, rows)
      priority[[length(priority) + 1L]] <<- f[order]
      modes[[length(modes) + 1L]] <<- list(title = title, kind = kind, target = target,
        reference = reference, method = method, total = nrow(table), rows = unname(tab[selected, , drop = FALSE]))
    }
    for (g in names(region$tables)) add_mode(region$tables[[g]], paste(g, "vs other regions"),
      "region", g, "other labelled regions", region$method)
    tests <- result$analysis[[name]]$comparisons$tests
    for (c in names(tests)) if (identical(tests[[c]]$status, "completed")) {
      t <- tests[[c]]
      tab <- t$table; tab$effect <- tab$logFC; tab$FDR <- tab$adj.P.Val
      add_mode(tab, paste("DE:", c), "test", t$contrast$numerator, t$contrast$denominator, t$method)
    }
    selected <- unique(as.integer(unlist(lapply(seq_len(min(features, max(c(lengths(priority), 0)))),
      function(i) vapply(priority, function(p) if (length(p) >= i) p[i] else NA_integer_, integer(1))))))
    selected <- head(unique(c(selected[!is.na(selected)], match(result$analysis[[name]]$pca$features, rownames(E)))), features)
    selected <- selected[!is.na(selected)]
    profiles <- lapply(selected, function(i) list(feature = i - 1L,
      values = vector(signif(as.numeric(E[i, ids]), 6))))
    native <- result$native_analysis$modalities[[name]]$samples
    graph <- lapply(names(native), function(s) {
      z <- native[[s]]
      if (is.null(z$graph_pca) || ncol(z$graph_pca$scores) < 2L) return(NULL)
      list(sample = s, coordinates = unname(z$graph_pca$scores[, 1:2, drop = FALSE]),
        pixels = vector(z$pixels), groups = vector(if (is.null(field)) rep(s, length(z$pixels)) else
          as.character(metadata[[field]][match(z$pixels, colnames(result$object))])))
    })
    payload$modalities[[length(payload$modalities) + 1L]] <- list(name = name,
      features = vector(rownames(E)), labels = vector(labels), modes = modes, profiles = profiles,
      region_names = vector(colnames(region$means)), region_counts = vector(region$groups$observations),
      region_means = if (is.null(region)) NULL else unname(signif(region$means, 6)),
      region_detected = if (is.null(region)) NULL else unname(signif(region$detected, 5)),
      graph_pca = Filter(Negate(is.null), graph))
  }
  payload
}

.wf_explorer_html <- function(result, points, features, rows) {
  data <- .wf_explorer_payload(result, points, features, rows)
  # Base64 keeps arbitrary feature/region names out of the HTML parser entirely.
  # Preserve small P/FDR values and cutoff decisions; only expression previews
  # are rounded explicitly in the payload builder.
  json <- jsonlite::toJSON(data, auto_unbox = TRUE, matrix = "rowmajor", digits = NA, na = "null", null = "null")
  encoded <- gsub("[\r\n]", "", jsonlite::base64_enc(charToRaw(enc2utf8(json))))
  paste0('<div id="spamtp-explorer" class="explorer">',
    '<div class="explorer-controls"><label>Modality<select id="ex-modality"></select></label>',
    '<label>Region / DE comparison<select id="ex-mode"></select></label>',
    '<label>Search features<input id="ex-search" type="search" placeholder="Feature ID or annotation"></label>',
    '<label>Direction<select id="ex-direction"><option value="all">Both</option><option value="up">Higher</option><option value="down">Lower</option></select></label>',
    '<label>Minimum |effect|<input id="ex-effect" type="number" min="0" step="0.1" value="0"></label>',
    '<label>Minimum directional AUC<input id="ex-auc" type="number" min="0.5" max="1" step="0.05" value="0.5"></label>',
    '<label>Minimum direction stability<input id="ex-stability" type="number" min="0" max="1" step="0.1" value="0"></label>',
    '<label>FDR cutoff<input id="ex-fdr" type="number" min="0" max="1" step="0.01" value="0.05"></label></div>',
    '<p id="ex-note" class="method-note" aria-live="polite"></p>',
    '<div class="explorer-plots"><div><h3 id="ex-effect-title">Region effects</h3><canvas id="ex-effects" width="640" height="360" aria-label="Effect plot; click a feature to inspect"></canvas></div>',
    '<div><div class="spatial-controls"><label>Sample<select id="ex-sample"></select></label>',
    '<label>Spatial colour<select id="ex-colour"><option value="expression">Selected feature</option><option value="region">Region membership</option></select></label>',
    '<button id="ex-reset" type="button">Reset view</button></div>',
    '<canvas id="ex-spatial" width="640" height="360" aria-label="Spatial feature map; scroll to zoom and drag to pan"></canvas></div></div>',
    '<h3 id="ex-feature-title">Selected feature</h3><p id="ex-profile-note" class="method-note"></p>',
    '<div class="explorer-plots"><div><h4>Mean expression across regions</h4><canvas id="ex-means" width="640" height="300" aria-label="Region mean expression"></canvas></div>',
    '<div><h4>Expression by region</h4><canvas id="ex-distribution" width="640" height="300" aria-label="Distribution of selected feature in the spatial preview"></canvas></div></div>',
    '<div class="explorer-actions"><button id="ex-download" type="button">Download filtered table</button><span id="ex-count"></span>',
    '<button id="ex-prev" type="button">Previous</button><button id="ex-next" type="button">Next</button></div>',
    '<div class="table-scroll"><table id="ex-table"><thead><tr><th><button data-sort="feature">Feature</button></th>',
    '<th><button data-sort="effect">Effect</button></th><th>Region mean</th><th>Reference mean</th>',
    '<th>Region non-zero fraction</th><th>Reference fraction</th><th><button data-sort="auc">Mean AUC</button></th>',
    '<th><button data-sort="cohen">Cohen effect</button></th><th>Direction stability</th><th><button data-sort="fdr">FDR</button></th></tr></thead><tbody></tbody></table></div>',
    '<div id="ex-tooltip" role="status" class="plot-tooltip" hidden></div>',
    '<noscript><p>Interactive browsing requires JavaScript. All static figures and complete CSV tables below remain available.</p></noscript>',
    '</div><script id="spamtp-explorer-data" type="application/json">"', encoded, '"</script>')
}

.wf_native_plots <- function(result, max_points) {
  plots <- list(); metadata <- as.data.frame(SummarizedExperiment::colData(result$object))
  region <- result$region_analysis$field
  xy <- SpatialExperiment::spatialCoords(result$object)
  spatial <- ncol(xy) >= 2L && all(is.finite(xy))
  for (name in names(result$analysis)) {
    assay <- result$settings$assay_names[[name]]
    features <- head(result$analysis[[name]]$pca$features, 3)
    if (!length(features)) next
    if (spatial) {
      ids <- .wf_preview_indices(metadata, region, max_points)
      small <- .wf_modality_spatial(result$object, assay, features, ids)
      panels <- plotSpatialFeature(small, features = features, assayName = "workflow",
        pointSize = 0.8, combine = FALSE)
      for (i in seq_along(panels)) plots[[paste0("native_map_", name, "_", i)]] <- panels[[i]] +
        ggplot2::labs(title = paste(name, features[i]), subtitle = "SpaMTP::plotSpatialFeature; recorded observation preview")
    }
    if (!is.null(region)) {
      counts <- sort(table(metadata[[region]]), decreasing = TRUE)
      groups <- head(names(counts)[counts >= 3 & nzchar(trimws(names(counts)))], 8)
      ids <- which(metadata[[region]] %in% groups)
      if (length(groups) >= 2 && length(ids) >= 4) {
        small <- result$object[, ids, drop = FALSE]
        plots[[paste0("native_violin_", name)]] <- mzViolinPlot(small, group.by = region,
          mzs = features, assay = assay, slot = "workflow", show.points = FALSE) +
          ggplot2::labs(title = paste(name, "regional feature distributions"),
            subtitle = "SpaMTP::mzViolinPlot; up to eight largest labelled groups") +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      }
    }
    samples <- result$native_analysis$modalities[[name]]$samples
    for (s in names(samples)) {
      a <- samples[[s]]
      if (!is.null(a$graph_pca) && ncol(a$graph_pca$scores) >= 2) {
        z <- a$graph_pca$scores
        frame <- data.frame(x = z[, 1], y = z[, 2], group = if (is.null(region)) s else
          as.character(metadata[[region]][match(a$pixels, colnames(result$object))]))
        plots[[paste0("native_graph_", name, "_", s)]] <- ggplot2::ggplot(frame,
          ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group)) +
          ggplot2::geom_point(size = 0.8, alpha = 0.65) +
          ggplot2::labs(title = paste(name, s, "spatial graph PCA"),
            subtitle = paste("SpaMTP::runSpatialGraphPCA;", nrow(z), "observations,", length(a$features), "features"))
      }
    }
    p <- result$analysis[[name]]$pathways
    if (!is.null(p$assay) && !is.null(region)) {
      E <- .assayData(result$object, p$assay, "pathwayScores")
      summary <- .wf_region_summary(E, metadata[[region]], 1)
      variability <- apply(summary$means, 1, stats::sd)
      selected <- head(order(variability, decreasing = TRUE, na.last = NA), 12)
      if (length(selected)) {
        frame <- expand.grid(pathway = rownames(E)[selected], region = colnames(summary$means))
        frame$score <- as.vector(summary$means[selected, , drop = FALSE])
        plots[[paste0("native_pathway_heatmap_", name)]] <- ggplot2::ggplot(frame,
          ggplot2::aes(x = .data$region, y = .data$pathway, fill = .data$score)) +
          ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#386b9d", mid = "white", high = "#be5c53") +
          ggplot2::labs(title = paste(name, "regional pathway scores"), subtitle = "createPathwayObject scores; top 12 variable region profiles") +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      }
    }
  }
  lapply(plots, function(p) p + ggplot2::theme_minimal(base_size = 11))
}

.wf_report_asset <- function(file) {
  path <- system.file("report", file, package = "SpaMTP")
  if (!nzchar(path)) stop("Missing installed report asset: ", file, call. = FALSE)
  paste(readLines(path, warn = FALSE), collapse = "\n")
}
