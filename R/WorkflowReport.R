.wf_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  gsub("'", "&#39;", x, fixed = TRUE)
}

.wf_flatten <- function(x) {
  x <- as.data.frame(x)
  for (name in names(x)) {
    if (is.list(x[[name]])) x[[name]] <- vapply(x[[name]], function(v)
      as.character(jsonlite::toJSON(v, auto_unbox = TRUE, null = "null", na = "null")), character(1))
    if (is.factor(x[[name]])) x[[name]] <- as.character(x[[name]])
  }
  x
}

.wf_html_table <- function(x, limit = 20) {
  x <- head(.wf_flatten(x), limit)
  if (!nrow(x)) return("<p>No rows.</p>")
  for (name in names(x)) if (is.numeric(x[[name]])) x[[name]] <- signif(x[[name]], 5)
  rows <- apply(x, 1, function(row) paste0("<tr><td>",
    paste(.wf_escape(row), collapse = "</td><td>"), "</td></tr>"))
  paste0('<div class="table-scroll"><table><thead><tr><th>',
    paste(.wf_escape(names(x)), collapse = "</th><th>"),
    "</th></tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table></div>")
}

.wf_report_plots <- function(x, max_points) {
  plots <- list()
  qc <- x$qc
  plots$qc <- ggplot2::ggplot(qc, ggplot2::aes(x = .data$modality, y = .data$detected)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "#a8c9df") +
    ggplot2::labs(x = NULL, y = "Non-zero features", title = "Input measurement coverage")
  n <- ncol(x$object)
  # Evenly spaced deterministic preview; analyses use all retained observations.
  idx <- unique(as.integer(round(seq(1, n, length.out = min(n, max_points)))))
  meta <- as.data.frame(SummarizedExperiment::colData(x$object))
  colour <- if (!is.null(x$settings$group)) as.character(meta[[x$settings$group]]) else
    if ("workflow_cluster" %in% names(meta)) as.character(meta$workflow_cluster) else rep("observations", n)
  embedding_plot <- function(E, title) {
    if (ncol(E) < 2L) return(NULL)
    frame <- data.frame(axis1 = E[idx, 1], axis2 = E[idx, 2], group = colour[idx])
    ggplot2::ggplot(frame, ggplot2::aes(x = .data$axis1, y = .data$axis2, colour = .data$group)) +
      ggplot2::geom_point(size = 0.65, alpha = 0.6) +
      ggplot2::labs(x = "PC 1", y = "PC 2", colour = x$settings$group %||% "Cluster",
        title = title, subtitle = paste(length(idx), "of", n, "retained observations displayed"))
  }
  for (name in names(x$analysis)) {
    pca <- x$analysis[[name]]$pca
    if (identical(pca$status, "completed"))
      plots[[paste0("pca_", name)]] <- embedding_plot(pca$scores, paste(name, "PCA"))
  }
  if (!is.null(x$joint)) plots$joint <- embedding_plot(x$joint$display, "Joint multi-omics structure")
  xy <- SpatialExperiment::spatialCoords(x$object)
  if (ncol(xy) >= 2L && all(is.finite(xy[, 1:2, drop = FALSE]))) {
    frame <- data.frame(x = xy[idx, 1], y = xy[idx, 2], group = colour[idx],
      sample = if ("sample_id" %in% names(meta)) as.character(meta$sample_id[idx]) else "sample")
    plots$spatial <- ggplot2::ggplot(frame, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group)) +
      ggplot2::geom_point(size = 0.65, alpha = 0.7) + ggplot2::coord_equal() +
      ggplot2::facet_wrap(~sample, scales = "free") +
      ggplot2::labs(title = "Retained observations in reference coordinates",
        subtitle = x$settings$coordinate_units %||% "Coordinate units were not declared")
    # Free scales and fixed aspect cannot be combined in recent ggplot2.
    plots$spatial <- plots$spatial + ggplot2::facet_wrap(~sample)
  }
  if (length(x$alignment)) {
    target <- .nativeCoordinates(x$object)
    for (name in names(x$alignment)) {
      a <- x$alignment[[name]]
      take <- function(frame) frame[unique(as.integer(round(seq(1, nrow(frame), length.out = min(nrow(frame), max_points))))), ]
      before <- take(a$before); after <- take(a$after); ref <- take(target)
      before$stage <- "Before"; after$stage <- "After"
      before$modality <- after$modality <- name
      ref$modality <- x$settings$reference
      r1 <- ref; r1$stage <- "Before"; r2 <- ref; r2$stage <- "After"
      fields <- c("x", "y", "modality", "stage")
      frame <- rbind(before[, fields], after[, fields], r1[, fields], r2[, fields])
      plots[[paste0("alignment_", name)]] <- ggplot2::ggplot(frame,
        ggplot2::aes(x = .data$x, y = .data$y, colour = .data$modality)) +
        ggplot2::geom_point(size = 0.65, alpha = 0.5) + ggplot2::coord_equal() +
        ggplot2::facet_wrap(~stage) + ggplot2::labs(title = paste(name, "registration overlay"),
          subtitle = "Geometric overlap alone does not establish biological correspondence")
    }
  }
  for (name in names(x$analysis)) {
    p <- x$analysis[[name]]$pathways
    if (!is.null(p)) {
      frame <- p$coverage[p$coverage$eligible, , drop = FALSE]
      if (!nrow(frame)) next
      plots[[paste0("coverage_", name)]] <- ggplot2::ggplot(frame,
        ggplot2::aes(x = .data$coverage_fraction)) + ggplot2::geom_histogram(bins = 30, fill = "#418a87") +
        ggplot2::labs(x = "Measured / database members", y = "Eligible pathways", title = paste(name, "pathway coverage"))
    }
  }
  lapply(plots, function(p) p + ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom"))
}

#' Render a portable, self-contained multi-omics report
#'
#' The HTML embeds PNG figures and styles and needs no network, Pandoc or R
#' session to view. Full tables, mapping weights, processing settings and the
#' analysed object are also saved for reproducibility. Existing files are never
#' overwritten. A report is published to output_dir only after rendering succeeds.
#' @param result A spamtp_workflow result from runSpaMTPWorkflow.
#' @param output_dir New or empty output directory.
#' @param max_points Maximum observations displayed per scatter plot. Analyses
#'   and exported tables are not subsampled.
#' @param interactive Embed an offline region/DE explorer using local JavaScript.
#' @param preview_features Maximum feature profiles embedded in the explorer.
#' @param preview_points Maximum spatial/distribution observations embedded.
#' @param table_rows Maximum preloaded ranked rows per region or contrast.
#'   Full tables are exported independently of this preview limit.
#' @return Absolute path to report.html, invisibly.
#' @export
#' @examples
#' utils::str(formals(renderSpaMTPReport))
renderSpaMTPReport <- function(result, output_dir, max_points = 5000,
    interactive = TRUE, preview_features = 180, preview_points = min(max_points, 2000), table_rows = 1000) {
  if (!inherits(result, "spamtp_workflow") || !identical(result$schema_version, 1L))
    stop("result must be a supported spamtp_workflow result.", call. = FALSE)
  .wf_number(max_points, "max_points", 3, TRUE)
  .wf_number(preview_features, "preview_features", 1, TRUE)
  .wf_number(preview_points, "preview_points", 3, TRUE)
  .wf_number(table_rows, "table_rows", 1, TRUE)
  if (!is.logical(interactive) || length(interactive) != 1 || is.na(interactive))
    stop("interactive must be TRUE or FALSE.", call. = FALSE)
  if (!is.character(output_dir) || length(output_dir) != 1L || !nzchar(output_dir))
    stop("output_dir must be a directory path.", call. = FALSE)
  if (file.exists(output_dir) && (!dir.exists(output_dir) ||
      length(list.files(output_dir, all.files = TRUE, no.. = TRUE))))
    stop("output_dir must be new or empty; existing reports are preserved.", call. = FALSE)
  dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
  output_dir <- file.path(normalizePath(dirname(output_dir), mustWork = TRUE), basename(output_dir))
  staging <- tempfile(".spamtp-report-", tmpdir = dirname(output_dir))
  dir.create(staging); on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  dir.create(file.path(staging, "tables")); dir.create(file.path(staging, "figures"))
  table_count <- 0L
  table_index <- data.frame(label = character(), file = character(), rows = integer())
  table <- function(data, label, preview = data) {
    data <- as.data.frame(data)
    table_count <<- table_count + 1L
    path <- sprintf("tables/%03d.csv", table_count)
    utils::write.csv(.wf_flatten(data), file.path(staging, path), row.names = FALSE, na = "")
    table_index[nrow(table_index) + 1L, ] <<- list(label, path, nrow(data))
    paste0("<h3>", .wf_escape(label), "</h3><p><a href=\"", path, "\">Download complete CSV (",
      nrow(data), " rows)</a></p>", .wf_html_table(preview))
  }
  basic <- .wf_report_plots(result, max_points)
  native_plots <- .wf_native_plots(result, max_points)
  if (!is.null(result$structure)) {
    basic <- basic[!grepl("^(pca_|joint|spatial)", names(basic))]
    native_plots <- native_plots[grepl("^native_pathway", names(native_plots))]
  }
  plot_html <- list(); plots <- c(basic, native_plots, .wf_analysis_plots(result, max_points))
  for (i in seq_along(plots)) {
    path <- file.path(staging, "figures", sprintf("%03d.png", i))
    size <- attr(plots[[i]], "report_size") %||% c(9, 5.5)
    caption <- attr(plots[[i]], "report_caption") %||% ""
    ggplot2::ggsave(path, plot = plots[[i]], width = size[1], height = size[2], dpi = 150, bg = "white")
    pdf <- sub("[.]png$", ".pdf", path)
    ggplot2::ggsave(pdf, plot = plots[[i]], width = size[1], height = size[2], bg = "white", device = grDevices::cairo_pdf)
    encoded <- jsonlite::base64_enc(readBin(path, "raw", n = file.info(path)$size))
    plot_html[[names(plots)[i]]] <- paste0('<figure><img width="', round(size[1] * 150),
      '" height="', round(size[2] * 150), '" alt="',
      .wf_escape(names(plots)[i]), '" src="data:image/png;base64,', encoded, '"><figcaption>',
      .wf_escape(caption), ' <a href="figures/', basename(pdf), '">PDF</a></figcaption></figure>')
  }
  plot_group <- function(prefix) paste(unlist(plot_html[startsWith(names(plot_html), prefix)]), collapse = "")
  widgets <- if (interactive) .wf_native_widgets(result, staging, min(max_points, 3000)) else
    list(html = list(associations = character(), pathways = character()), status = data.frame())
  section <- function(id, title, content) paste0('<section id="', id, '"><h2>', title, '</h2>', content, '</section>')
  specs <- result$settings$modalities
  inventory <- do.call(rbind, lapply(names(specs), function(n) {
    s <- specs[[n]]; e <- .experimentForAssay(result$object, result$settings$assay_names[[n]])
    data.frame(modality = n, type = s$type, features = nrow(e), retained_observations = ncol(e),
      input_layer = s$layer, normalization = s$normalization, species = s$species %||% "undeclared")
  }))
  overview <- paste0('<p class="lead">', length(specs), ' modalities &middot; ', ncol(result$object),
    ' retained observations &middot; ', .wf_escape(result$settings$observation_unit), '</p>',
    '<p>Follow the measurement and alignment evidence before interpreting shared structure. ',
    'The report distinguishes computed, provided, skipped and unrequested stages.</p>',
    table(inventory, "Input and processing choices"),
    '<details><summary>Processing stages and methods</summary>', table(result$stages, "Workflow status"), '</details>')
  contract <- data.frame(
    analytical_question = c("What is measured and retained?", "What changes after adding space or another modality?",
      "Which features characterize a region?", "What is supported across biological replicates?",
      "Does cross-omic association survive regional composition adjustment?", "Which pathway members were actually measured?"),
    analysis = c("Declared layers, QC, correspondence weights and source provenance",
      "Matched-feature PCA / Graph PCA / multiOmicIntegration; K sensitivity and external-reference ARI",
      "findAllDEMs marker effects, pairwise AUC, native heatmaps and spatial-block sensitivity",
      "findAllDEMs replicate mode: explicit paired or independent design, full test family and effect intervals",
      "findCorrelatedFeatures: raw and covariate-residual correlations; paired spatial layers",
      "One identity assay/index for complete-rank enrichment, scores, GESECA and network views"),
    inference_unit = c(result$settings$observation_unit, "Exploratory tissue structure; reference labels used only for evaluation",
      "Observed region distributions; no pixel-based P values",
      result$settings$replicate %||% "No biological replicate column supplied",
      "Descriptive matched-observation association", "Competitive feature-set test; measured panel limits apply"))
  overview <- paste0(overview, table(contract, "Analysis contract"))
  if (length(result$warnings)) overview <- paste0(overview,
    '<details open><summary>Warnings recorded during analysis</summary><ul>',
    paste0('<li>', .wf_escape(result$warnings), '</li>', collapse = ""), '</ul></details>')
  qcsummary <- do.call(rbind, lapply(split(result$qc, result$qc$modality), function(q)
    data.frame(modality = q$modality[1], input_observations = nrow(q), qc_pass = sum(q$pass),
      excluded = sum(!q$pass), median_detected = stats::median(q$detected))))
  qcbody <- paste0(table(qcsummary, "QC summary"), plot_group("qc"), table(result$qc, "Per-observation QC"),
    table(result$retained, "Joint inclusion audit"),
    '<p>QC thresholds apply to non-zero features. Constant features remain in the measured pathway universe. ',
    'The workflow assay is the selected layer after the stated normalization; original layers are retained. ',
    'Library scaling does not correct batch effects.</p>')
  alignmentbody <- if (length(result$mapping)) paste0(plot_group("alignment_"),
    paste(vapply(names(result$mapping), function(n) {
      m <- result$mapping[[n]]; used <- lengths(m$matches) > 0
      paste0('<p>', .wf_escape(n), ': ', sum(used), '/', length(used),
        ' reference observations mapped; ', length(unique(unlist(m$matches))),
        ' unique source observations contributed. Reused source measurements are not independent replicates.</p>',
        table(m$observations, paste(n, "mapping coverage")))
    }, character(1)), collapse = "")) else
      '<p>Pairing was supplied by the input alternative experiments. This run did not estimate or validate a registration. Provenance from the input is retained in workflow.rds.</p>'
  if (ncol(SpatialExperiment::spatialCoords(result$object)) < 2L)
    alignmentbody <- paste0(alignmentbody, '<p>No usable spatial coordinates were supplied; no spatial map is drawn.</p>')
  singlebody <- paste0(plot_group("pca_"),
    '<p>PCA uses variable features of each normalized modality. Feature selection and variance explained are retained in workflow.rds.</p>')
  jointbody <- if (!is.null(result$joint)) paste0('<p>', .wf_escape(result$joint$method),
    '</p>', plot_group("joint"), plot_group("spatial"),
    '<p>Structure and optional clusters describe this input. They are not validated cell types, mechanisms or evidence of cross-sample reproducibility.</p>') else
      paste0('<p>Joint embedding was not estimable; see workflow status.</p>', plot_group("spatial"))
  associationbody <- paste0('<p>These are exploratory correlations across retained observations. Spatial autocorrelation, ',
    'shared group effects and source-pixel reuse can drive them. No pixel-level P values are calculated. ',
    'Only the configured most variable features enter this screen.</p>',
    paste(vapply(names(result$associations), function(n)
      table(result$associations[[n]], n), character(1)), collapse = ""))
  if (!is.null(result$structure)) {
    singlebody <- paste0('<p>', .wf_escape(result$structure$interpretation), '</p>',
      plot_group("structure_embedding"), plot_group("structure_spatial"))
    jointbody <- paste0('<p>The primary representation and K are specified before evaluation. ',
      'Reference annotation is used only for evaluation. Inspect performance across K and disagreement between modalities.</p>',
      plot_group("structure_evaluation"), table(result$structure$metrics, "Structure comparison across K"),
      table(result$structure$blocks, "Spatial blocks used for conditional clustering sensitivity"))
    associationbody <- paste0('<p>Regional marker anchors are queried with findCorrelatedFeatures(). ',
      'Raw correlations and correlations after regressing the declared region/sample/replicate fields ',
      'answer whether association extends beyond a shared regional shift. Residual correlation does not establish mechanism.</p>',
      plot_group("association_adjustment_"), plot_group("association_pair_"),
      '<details><summary>Complete screened pairs and conditional analysis settings</summary>',
      paste(vapply(names(result$associations), function(n) paste0(table(result$associations[[n]], n),
        '<pre>', .wf_escape(paste(capture.output(str(attr(result$associations[[n]], "association"), max.level = 1)), collapse = "\n")), '</pre>'), character(1)), collapse = ""), '</details>')
  }
  regionbody <- if (interactive && !is.null(result$region_analysis))
    .wf_explorer_html(result, preview_points, preview_features, table_rows) else
      '<p class="method-note">Interactive browsing is disabled or region summaries are absent. Use analyzeSpaMTPRegions() to add summaries to an earlier workflow result.</p>'
  regiontables <- character()
  for (n in names(result$region_analysis$modalities)) {
    r <- result$region_analysis$modalities[[n]]
    regiontables <- c(regiontables, table(r$groups, paste(n, "region eligibility")))
    if (length(r$tables)) {
      combined <- do.call(rbind, lapply(names(r$tables), function(g) data.frame(region = g, r$tables[[g]])))
      regiontables <- c(regiontables, table(combined, paste(n, "all region effects")))
    }
    if (r$unlabelled) regiontables <- c(regiontables, paste0('<p class="method-note">', r$unlabelled,
      ' observations have blank/missing region labels and are excluded from region effects.</p>'))
  }
  regionbody <- paste0(regionbody, '<details><summary>Complete regional tables and eligibility</summary>',
    paste(regiontables, collapse = ""), '</details>')
  regionbody <- paste0('<p>Characterize each region with pairwise effect sizes and detection, then inspect ',
    'the spatial distributions of the same ranked features. Region characterization is conditional on the supplied or learned partition.</p>',
    plot_group("markers_heatmap_"), plot_group("markers_detection_"), plot_group("markers_stability_"),
    regionbody, plot_group("markers_maps_"))
  nativebody <- paste0('<p class="method-note">The following outputs are computed by SpaMTP native functions. ',
    'Spatial statistics and graph PCA use the recorded bounded feature/point screen, separately by sample. ',
    'These exploratory spatial statistics do not provide biological-replicate DE evidence.</p>',
    plot_group("native_map_"), plot_group("native_violin_"), plot_group("native_graph_"))
  for (n in names(result$native_analysis$modalities)) {
    samples <- result$native_analysis$modalities[[n]]$samples
    for (s in names(samples)) {
      a <- samples[[s]]
      if (!is.null(a$moran)) nativebody <- paste0(nativebody,
        table(a$moran, paste(n, s, "Moran's I (native function)"),
          a$moran[order(a$moran$FDR), , drop = FALSE]))
      if (identical(a$status, "skipped")) nativebody <- paste0(nativebody,
        '<p class="method-note">', .wf_escape(paste(n, s, a$reason)), '</p>')
    }
  }
  capabilities <- data.frame(
    capability = c("Region effects and linked DE browsing", "Spatial feature maps", "Regional feature distributions",
      "Spatial autocorrelation", "Spatial graph PCA", "Joint multi-omics", "Mass annotation", "Pathway scores and coverage", "Measured-universe enrichment"),
    function_name = c("analyzeSpaMTPRegions", "plotSpatialFeature", "mzViolinPlot", "findSpatiallyVariableMetabolites",
      "runSpatialGraphPCA", "multiOmicIntegration", "annotateSM", "createPathwayObject / buildPathwayIndex", "fishersPathwayAnalysis"),
    status = c(result$region_analysis$status %||% "not available",
      if (any(startsWith(names(plots), "native_map_"))) "shown" else "no coordinates",
      if (any(startsWith(names(plots), "native_violin_"))) "shown" else "no eligible groups",
      if (any(vapply(result$native_analysis$modalities, function(m) any(vapply(m$samples, function(s) !is.null(s$moran), logical(1))), logical(1)))) "computed" else "not run",
      if (any(startsWith(names(plots), "native_graph_"))) "computed" else "not run",
      if (is.null(result$joint)) "not applicable" else "computed",
      if (any(vapply(result$analysis, function(a) !is.null(a$annotation), logical(1)))) "computed" else "not requested",
      if (any(vapply(result$analysis, function(a) !is.null(a$pathways$assay), logical(1)))) "computed" else "no configured index",
      if (any(vapply(result$analysis, function(a) !is.null(a$pathways$enrichment), logical(1)))) "computed" else "no foreground"))
  cards <- vapply(seq_len(nrow(capabilities)), function(i) paste0('<div class="capability"><span class="status">',
    .wf_escape(capabilities$status[i]), '</span><strong>', .wf_escape(capabilities$capability[i]),
    '</strong><p class="method-note function-name">', .wf_escape(capabilities$function_name[i]), '</p></div>'), character(1))
  if (is.null(result$structure)) overview <- paste0(overview, '<h3>Native functionality in this report</h3><div class="capability-grid">', paste(cards, collapse = ""), '</div>')
  comparisons <- character()
  for (n in names(result$analysis)) {
    comp <- result$analysis[[n]]$comparisons
    if (!is.null(comp$group_means)) {
      means <- data.frame(feature = rownames(comp$group_means), comp$group_means, check.names = FALSE)
      comparisons <- c(comparisons, table(means, paste(n, "descriptive group means")))
    }
    for (name in names(comp$tests)) {
      t <- comp$tests[[name]]
      comparisons <- c(comparisons, paste0('<h3>', .wf_escape(paste(n, name)), '</h3><p>',
        .wf_escape(t$method %||% t$reason), '</p>'), table(t$units, "Biological replicate design"))
      if (!is.null(t$table)) comparisons <- c(comparisons, table(t$table, "Replicate-level contrast",
        t$table[order(t$table$adj.P.Val), , drop = FALSE]))
    }
  }
  contrastbody <- paste0('<p>Without a configured biological-replicate contrast, group means are descriptive. ',
    'Inferential results compare equally weighted replicate means of processed values using limma; ',
    'the logFC column is a difference on the workflow scale, not necessarily a raw abundance fold change. ',
    'BH correction is within each modality and contrast. Paired tests use complete pairs.</p>',
    '<details><summary>Group means, replicate design and full DE results</summary>', paste(comparisons, collapse = ""), '</details>')
  pathwaybody <- paste0('<p>Mass matches are putative annotations. Gene pathway scores describe measured members, not biochemical activity. Missing identities are never counted as database IDs.</p>', plot_group("native_pathway_heatmap_"))
  contrastbody <- paste0(plot_group("replicate_effects_"), plot_group("replicate_units_"), contrastbody)
  for (n in names(result$analysis)) {
    a <- result$analysis[[n]]
    if (is.null(a$pathways)) pathwaybody <- paste0(pathwaybody,
      '<p>', .wf_escape(n), ': pathway analysis was not configured. A compatible identity and pathway index is required for ',
      .wf_escape(result$settings$modalities[[n]]$species %||% 'the declared species'), '.</p>')
    if (!is.null(a$annotation)) pathwaybody <- paste0(pathwaybody,
      table(a$annotation, paste(n, "mass annotations")))
    if (!is.null(a$pathways)) {
      p <- a$pathways
      if (identical(p$status, 'skipped')) pathwaybody <- paste0(pathwaybody,
        '<p>', .wf_escape(p$reason), '</p>')
      columns <- intersect(c("pathwayRampId", "pathwayName", "database_size", "database_raw_size",
        "measured_size", "coverage_fraction", "used_size", "eligible"), names(p$coverage))
      pathwaybody <- paste0(pathwaybody, plot_group(paste0("coverage_", n)),
        table(p$coverage, paste(n, "complete pathway coverage and membership"), p$coverage[p$coverage$eligible, columns, drop = FALSE]),
        table(p$mapping, paste(n, "gene identity audit")))
      if (!is.null(p$enrichment)) pathwaybody <- paste0(pathwaybody,
        '<p>Enrichment uses the supplied foreground and the measured modality as universe, including eligible zero-overlap pathways in the correction family.</p>',
        table(p$enrichment, paste(n, "foreground enrichment")))
    }
  }
  provenance <- paste0('<p>Core analysis created ', .wf_escape(result$created), ' &middot; SpaMTP ', .wf_escape(result$package_version),
    ' &middot; seed ', result$settings$seed, '</p><p>Region/native extension: SpaMTP ',
    .wf_escape(result$region_analysis$package_version %||% "unavailable"),
    ' &middot; report renderer: SpaMTP ', as.character(utils::packageVersion("SpaMTP")),
    '</p><p><a href="workflow.rds">Complete workflow RDS</a> &middot; ',
    '<a href="sessionInfo.txt">R session</a> &middot; <a href="manifest.csv">File checksums</a></p>',
    '<details><summary>Processing configuration</summary><pre>',
    .wf_escape(paste(capture.output(str(result$settings, max.level = 4)), collapse = "\n")),
    '</pre></details><p>The RDS retains source metadata, transformations, mapping weights, pathway index provenance, ',
    'complete mapping conflicts, analysis matrices and method settings. ',
    'Data are not sent to external services by rendering. Remote input acquisition only occurs when explicitly configured.</p>')
  associationbody <- paste0(associationbody, paste(widgets$html$associations, collapse = ""))
  pathwaybody <- paste0(pathwaybody, plot_group("pathway_regional_"), plot_group("pathway_geseca_"),
    paste(widgets$html$pathways, collapse = ""))
  for (name in names(result$analysis)) {
    p <- result$analysis[[name]]$pathways
    if (!is.null(p$regional)) {
      equivalence <- .wf_pathway_equivalence(p$regional)
      if (nrow(equivalence)) {
        result$analysis[[name]]$pathways$member_equivalence <- equivalence
        pathwaybody <- paste0(pathwaybody, '<p>', .wf_escape(name), ': ',
          length(unique(p$regional$pathwayRampId)), ' tested pathway labels represent ', nrow(equivalence),
          ' distinct patterns of measured members across regions. Labels sharing the same measured members do not provide independent evidence. ',
          'The complete pathway family remains in BH correction.</p>',
          '<details><summary>Measured-member equivalence audit</summary>',
          table(equivalence, paste(name, "measured pathway redundancy")), '</details>')
      }
      pathwaybody <- paste0(pathwaybody,
        table(p$regional, paste(name, "complete regional pathway family")))
    }
    if (!is.null(p$geseca)) pathwaybody <- paste0(pathwaybody,
      table(p$geseca, paste(name, "complete GESECA results")))
    if (!is.null(p$compound_mapping)) pathwaybody <- paste0(pathwaybody,
      table(p$compound_mapping$inputs, paste(name, "compound candidates before ambiguity exclusion")),
      '<p>', length(p$compound_mapping$excluded_conflicts), ' ambiguous mass features excluded from pathway identities.</p>')
  }
  provenance <- paste0(provenance, table(widgets$status, "Native interactive views"))
  result$native_widget_status <- widgets$status
  if (!is.null(result$structure) && any(vapply(result$native_analysis$modalities,
      function(m) any(vapply(m$samples, function(s) !is.null(s$moran), logical(1))), logical(1))))
    jointbody <- paste0(jointbody, '<details><summary>Optional native spatial autocorrelation screen</summary>',
      nativebody, '</details>')
  content <- paste0(section("overview", "Study and analysis design", overview),
    section("qc", "Measurement quality and preparation", qcbody),
    section("alignment", "Registration and correspondence", alignmentbody),
    section("single", "Individual, spatial and joint representations", singlebody),
    section("joint", "Structure evaluation and sensitivity", jointbody),
    section("regions", "Region markers: effect, specificity and stability", regionbody),
    if (is.null(result$structure)) section("native", "Additional native spatial analyses", nativebody) else "",
    section("contrasts", "Biological-replicate contrasts", contrastbody),
    section("associations", "Cross-omic colocalization and regional composition", associationbody),
    section("pathways", "Identity, annotation and pathways", pathwaybody),
    section("methods", "Methods and reproducibility", provenance))
  # Explicit method notes and figure captions use smaller serif text; narrative
  # results and navigation retain the readable body face.
  style <- .wf_report_asset("report.css")
  html <- paste0('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>',
    .wf_escape(result$title), '</title><style>', style, '</style></head><body><header><h1>', .wf_escape(result$title),
    '</h1><p>Acquisition &rarr; QC &rarr; alignment &rarr; analysis &rarr; interpretation</p><nav>',
    paste0('<a href="#', c("overview", "qc", "alignment", "single", "joint", "regions", "contrasts", "associations", "pathways", "methods"), '">',
      c("Design", "QC / histology", "Correspondence", "Representations", "Evaluation", "Region markers / DE", "Replicate contrasts", "Colocalization", "Pathways", "Reproduce"), '</a>', collapse = ""),
    '</nav></header><main>', content, '</main>',
    if (interactive) paste0('<script type="text/javascript">', .wf_report_asset("explorer.js"), '</script>') else "", '</body></html>')
  writeLines(html, file.path(staging, "report.html"), useBytes = TRUE)
  utils::write.csv(table_index, file.path(staging, "table_index.csv"), row.names = FALSE)
  utils::write.csv(data.frame(name = names(plots), png = sprintf("figures/%03d.png", seq_along(plots)),
    pdf = sprintf("figures/%03d.pdf", seq_along(plots)), caption = vapply(plots,
      function(p) attr(p, "report_caption") %||% "", character(1))),
    file.path(staging, "figure_index.csv"), row.names = FALSE)
  result$report <- file.path(output_dir, "report.html")
  result$report_settings <- list(interactive = interactive, max_points = max_points,
    preview_points = preview_points, preview_features = preview_features, table_rows = table_rows,
    renderer_version = as.character(utils::packageVersion("SpaMTP")))
  saveRDS(result, file.path(staging, "workflow.rds"))
  writeLines(c("Core analysis session:", capture.output(result$session),
    "", "Report rendering session:", capture.output(utils::sessionInfo())), file.path(staging, "sessionInfo.txt"))
  files <- list.files(staging, recursive = TRUE)
  utils::write.csv(data.frame(file = files, bytes = file.info(file.path(staging, files))$size,
    md5 = unname(tools::md5sum(file.path(staging, files)))), file.path(staging, "manifest.csv"), row.names = FALSE)
  if (dir.exists(output_dir)) unlink(output_dir, recursive = FALSE)
  if (!file.rename(staging, output_dir)) stop("Could not publish rendered report directory.", call. = FALSE)
  invisible(file.path(output_dir, "report.html"))
}
