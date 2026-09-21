.wf_figure <- function(plot, width = 12, height = 7, caption = "") {
  attr(plot, "report_size") <- c(width, height)
  attr(plot, "report_caption") <- caption
  plot
}

.wf_plot_theme <- function() ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(size = 10, colour = "#526574"),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "bottom", legend.key.height = grid::unit(3, "mm"))

.wf_analysis_plots <- function(result, max_points) {
  plots <- list(); obj <- result$object; n <- ncol(obj)
  field <- result$region_analysis$field
  md <- as.data.frame(SummarizedExperiment::colData(obj))
  idx <- .wf_preview_indices(md, field, max_points)
  point_size <- if (length(idx) < 300L) 1.8 else if (length(idx) < 1500L) 1.2 else 0.85
  xy <- SpatialExperiment::spatialCoords(obj)
  spatial <- all(c("x", "y") %in% colnames(xy)) && all(is.finite(xy))
  reps <- result$structure$representations
  reference <- result$structure$config$reference %||% field
  colour <- if (is.null(reference)) as.character(obj$sample_id) else as.character(obj[[reference]])
  levels <- sort(unique(stats::na.omit(colour)))
  palette <- stats::setNames(grDevices::hcl.colors(length(levels), "Dark 3"), levels)
  images <- SpatialExperiment::imgData(obj)
  if (spatial && nrow(images)) {
    panels <- lapply(seq_len(nrow(images)), function(i) {
      sample <- images$sample_id[i]
      image <- .nativeImage(obj, images$image_id[i], sample)
      take <- idx[as.character(obj$sample_id[idx]) == sample]
      frame <- data.frame(x = xy[take, "x"], y = xy[take, "y"], region = colour[take])
      ggplot2::ggplot(frame, ggplot2::aes(.data$x, .data$y)) + .imageLayer(image, 0.85) +
        ggplot2::geom_point(ggplot2::aes(colour = .data$region), size = point_size, alpha = .75) +
        ggplot2::scale_colour_manual(values = palette) +
        ggplot2::scale_y_reverse() + ggplot2::coord_equal() + .wf_plot_theme() +
        ggplot2::labs(title = as.character(sample), x = NULL, y = NULL, colour = reference %||% "Sample")
    })
    plots$qc_histology <- .wf_figure(cowplot::plot_grid(plotlist = panels, ncol = min(2, length(panels))),
      caption = "Registered imgData raster and scaleFactor in the same coordinate frame as the observations. Histology was acquired from the recorded source archive; no image alignment was inferred from expression.")
  }
  if (length(reps)) {
    frame <- do.call(rbind, lapply(names(reps), function(name) {
      r <- reps[[name]]
      data.frame(x = r$display[idx, 1], y = r$display[idx, 2], colour = colour[idx],
        modality = r$modality, method = r$method)
    }))
    frame$modality <- factor(frame$modality, levels = c(names(result$analysis), "Joint"))
    frame$method <- factor(frame$method, levels = c("PCA", "Graph PCA"))
    plots$structure_embedding <- .wf_figure(ggplot2::ggplot(frame,
      ggplot2::aes(.data$x, .data$y, colour = .data$colour)) +
      ggplot2::geom_point(size = point_size, alpha = 0.8) +
      ggplot2::scale_colour_manual(values = palette) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))) +
      ggplot2::facet_grid(method ~ modality, scales = "free") +
      ggplot2::labs(title = "Does joint structure agree with either individual modality?",
        subtitle = "Each panel uses the full representation; the same annotation colours apply across panels",
        x = "Embedding 1", y = "Embedding 2", colour = reference %||% "Sample") + .wf_plot_theme(),
      height = 7.5, caption = paste("Display:", paste(unique(vapply(reps,
        function(r) r$display_method, character(1))), collapse = "; "),
        ". Clustering and evaluation use all retained dimensions, not these two axes. Embedding orientations and distances between panels are not comparable."))
    k <- result$settings$clusters
    if (!is.null(k) && spatial) {
      frame <- do.call(rbind, lapply(names(reps), function(name) {
        r <- reps[[name]]; labs <- result$structure$assignments[[name]][[as.character(k)]]
        if (is.null(labs)) return(NULL)
        data.frame(x = xy[idx, "x"], y = xy[idx, "y"], cluster = factor(labs[idx]),
          modality = r$modality, method = r$method, sample = obj$sample_id[idx])
      }))
      frame$modality <- factor(frame$modality, levels = c(names(result$analysis), "Joint"))
      frame$method <- factor(frame$method, levels = c("PCA", "Graph PCA"))
      plots$structure_spatial <- .wf_figure(ggplot2::ggplot(frame,
        ggplot2::aes(.data$x, .data$y, colour = .data$cluster)) +
        ggplot2::geom_point(size = point_size) + ggplot2::scale_y_reverse() + ggplot2::coord_equal() +
        ggplot2::scale_colour_manual(values = stats::setNames(grDevices::hcl.colors(k, "Dark 3"), seq_len(k))) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3))) +
        ggplot2::facet_grid(method + sample ~ modality) +
        ggplot2::labs(title = paste("Where do the partitions differ?  K =", k),
          subtitle = "Identical spatial coordinates and axis limits; cluster numbers are local to each panel",
          x = NULL, y = NULL, colour = "Local cluster") + .wf_plot_theme(), height = 8,
        caption = "Graph regularization encourages neighboring observations to have similar embeddings. Smooth spatial domains alone do not demonstrate better recovery of tissue anatomy.")
    }
    metrics <- result$structure$metrics
    if (nrow(metrics)) {
      selected <- c("approx_silhouette", "reference_ARI", "conditional_stability", "spatial_neighbor_agreement")
      frame <- do.call(rbind, lapply(selected, function(m) data.frame(
        representation = metrics$representation, k = metrics$k, metric = m, value = metrics[[m]])))
      frame <- frame[is.finite(frame$value), , drop = FALSE]
      if (nrow(frame)) plots$structure_evaluation <- .wf_figure(ggplot2::ggplot(frame,
        ggplot2::aes(.data$k, .data$value, colour = .data$representation, group = .data$representation)) +
        ggplot2::geom_line(linewidth = 0.65) + ggplot2::geom_point(size = 2) +
        ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 2) +
        ggplot2::scale_x_continuous(breaks = unique(metrics$k)) +
        ggplot2::labs(title = "Sensitivity to the requested number of clusters", x = "K", y = NULL, colour = "Representation") + .wf_plot_theme(),
        height = 8, caption = result$structure$interpretation)
    }
  }
  for (name in names(result$region_analysis$modalities)) {
    m <- result$region_analysis$modalities[[name]]$markers
    if (!identical(m$status, "completed")) next
    tab <- m$DEMs
    top <- do.call(rbind, lapply(split(tab, tab$cluster), function(t)
      head(t[order(-t$mean_auc, t$gene), , drop = FALSE], 3)))
    heat <- demsHeatmap(m, n = 3, only.pos = TRUE, logfc.threshold = 0,
      order.by = "mean_auc", fontsize_row = 8, fontsize_col = 10,
      color = grDevices::colorRampPalette(c("#3b5b92", "#f7f7f3", "#b73e45"))(101),
      annotation_colors = palette)
    plots[[paste0("markers_heatmap_", name)]] <- .wf_figure(cowplot::ggdraw(heat$gtable),
      width = 11, height = max(5, length(heat$selected_features) * 0.23 + 2),
      caption = paste(name, "\u2014 demsHeatmap(): top three positive markers per eligible region, ranked by mean pairwise AUC. Colour is row-scaled mean workflow expression. The displayed features and exact unscaled matrix are recorded in the marker results."))
    features <- unique(top$gene)
    means <- m$expression[features, , drop = FALSE]
    z <- t(scale(t(means))); z[!is.finite(z)] <- 0
    frame <- do.call(rbind, lapply(seq_len(ncol(means)), function(i) {
      t <- tab[tab$cluster == colnames(means)[i], ]
      data.frame(feature = features, region = colnames(means)[i], z = z[, i],
        detection = t$detected_region[match(features, t$gene)])
    }))
    frame$feature <- factor(frame$feature, levels = rev(features))
    plots[[paste0("markers_detection_", name)]] <- .wf_figure(ggplot2::ggplot(frame,
      ggplot2::aes(.data$region, .data$feature, size = .data$detection, colour = .data$z)) +
      ggplot2::geom_point() + ggplot2::scale_colour_gradient2(low = "#3b5b92", mid = "#f7f7f3", high = "#b73e45") +
      ggplot2::scale_size_area(max_size = 6, limits = c(0, 1)) +
      ggplot2::labs(title = paste(name, "\u2014 are markers enriched or merely widely detected?"),
        x = "Region", y = NULL, size = "Non-zero fraction", colour = "Mean z-score") + .wf_plot_theme() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1)),
      height = max(5, length(features) * 0.22 + 2))
    stable <- top[is.finite(top$cohen_min_block) & is.finite(top$cohen_max_block) & is.finite(top$mean_cohen), ]
    if (nrow(stable)) {
      stable$label <- factor(paste(stable$cluster, stable$gene, sep = " | "),
        levels = rev(unique(paste(stable$cluster, stable$gene, sep = " | "))))
      plots[[paste0("markers_stability_", name)]] <- .wf_figure(ggplot2::ggplot(stable,
        ggplot2::aes(.data$mean_cohen, .data$label, colour = .data$mean_auc)) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey75") +
        ggplot2::geom_segment(ggplot2::aes(x = .data$cohen_min_block, xend = .data$cohen_max_block,
          yend = .data$label), linewidth = 0.9) + ggplot2::geom_point(size = 2) +
        ggplot2::scale_colour_viridis_c(limits = c(0.5, 1), oob = scales::squish) +
        ggplot2::labs(title = paste(name, "\u2014 does a spatial subregion drive the marker?"),
          subtitle = "Line: range after leaving one spatial block out; point: full-tissue estimate",
          x = "Mean pairwise Cohen effect", y = NULL, colour = "Mean AUC") + .wf_plot_theme(),
        height = max(5, nrow(stable) * .22 + 2),
        caption = "These ranges describe sensitivity to coordinate-only block omission. They are not confidence intervals. Regions are fixed at their full-data labels; rare regions may prevent some omissions, as recorded in the audit.")
    }
    if (spatial) {
      chosen <- vapply(split(tab, tab$cluster), function(t) t$gene[order(-t$mean_auc, t$gene)[1]], character(1))
      chosen <- unique(chosen)
      panel <- if (nrow(images)) spatialMZPlot(obj[, idx], mzs = chosen,
        assay = result$settings$assay_names[[name]], slot = "workflow", images = unique(images$image_id),
        combine = FALSE, pt.size.factor = point_size, image.alpha = .6,
        min.cutoff = "q02", max.cutoff = "q98", keep.scale = "none") else
        imageMZPlot(obj[, idx], mzs = chosen, assay = result$settings$assay_names[[name]],
          slot = "workflow", combine = FALSE, dark.background = FALSE, size = point_size,
          cols = c("#eef0f1", "#714068"), min.cutoff = "q02", max.cutoff = "q98", scale = "none")
      panel <- lapply(panel, function(p) {
        if (!is.null(p$scales$get_scales("fill"))) p <- p + ggplot2::labs(fill = "Value")
        if (!is.null(p$scales$get_scales("colour"))) p <- p + ggplot2::labs(colour = "Value")
        p + ggplot2::theme(plot.title = ggplot2::element_text(size = 11),
          legend.title = ggplot2::element_text(size = 9),
          legend.text = ggplot2::element_text(size = 8), legend.position = "bottom")
      })
      plots[[paste0("markers_maps_", name)]] <- .wf_figure(
        cowplot::plot_grid(plotlist = panel, ncol = 3), height = ceiling(length(panel) / 3) * 3.3,
        caption = paste(name, if (nrow(images)) "\u2014 spatialMZPlot():" else "\u2014 imageMZPlot():",
          "one highest-AUC feature per eligible region. Colours are independently clipped at the 2nd/98th percentiles; numeric intensities are not comparable across features. All panels retain the same tissue coordinates."))
    }
  }
  for (name in names(result$associations)) {
    t <- result$associations[[name]]; finite <- is.finite(t$correlation) & is.finite(t$correlation_raw)
    d <- t[finite, , drop = FALSE]
    if (!nrow(d)) next
    plots[[paste0("association_adjustment_", name)]] <- .wf_figure(ggplot2::ggplot(d,
      ggplot2::aes(.data$correlation_raw, .data$correlation)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey80") + ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70") +
      ggplot2::geom_point(size = 0.65, alpha = 0.3, colour = "#256e87") +
      ggplot2::coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
      ggplot2::labs(title = paste(name, "\u2014 how much correlation is explained by region and sample?"),
        x = "Raw Pearson correlation", y = "Residual Pearson correlation") + .wf_plot_theme(),
      caption = paste(attr(t, "association")$interpretation, "Covariates:", paste(attr(t, "association")$covariates, collapse = ", ")))
    best <- d[1, ]; E1 <- .assayData(obj, result$settings$assay_names[[best$modality1]], "workflow")
    E2 <- .assayData(obj, result$settings$assay_names[[best$modality2]], "workflow")
    ids <- match(attr(t, "association")$observations, colnames(obj))
    points <- data.frame(x = as.numeric(E1[best$feature1, ids]), y = as.numeric(E2[best$feature2, ids]), region = colour[ids])
    plots[[paste0("association_pair_", name)]] <- .wf_figure(ggplot2::ggplot(points,
      ggplot2::aes(.data$x, .data$y, colour = .data$region)) +
      ggplot2::geom_point(size = point_size, alpha = 0.5) + ggplot2::scale_colour_manual(values = palette) +
      ggplot2::labs(title = "Inspect the highest-ranked conditional association",
        subtitle = sprintf("Raw r = %.3f; residual r = %.3f; selected from %d screened pairs", best$correlation_raw, best$correlation, nrow(t)),
        x = paste(best$modality1, best$feature1), y = paste(best$modality2, best$feature2), colour = reference %||% "Sample") + .wf_plot_theme(),
      caption = "This pair is selected on these data and requires independent validation. The points show the original workflow values, with the same annotation colours used in the structure panels.")
  }
  for (name in names(result$analysis)) for (contrast in names(result$analysis[[name]]$comparisons$tests)) {
    t <- result$analysis[[name]]$comparisons$tests[[contrast]]
    if (!identical(t$status, "completed") || is.null(t$expression)) next
    top <- head(t$table[order(t$table$adj.P.Val, t$table$feature), ], 12)
    plots[[paste0("replicate_effects_", name, "_", contrast)]] <- .wf_figure(ggplot2::ggplot(top,
      ggplot2::aes(.data$logFC, stats::reorder(.data$feature, .data$logFC), colour = .data$adj.P.Val < .05)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
      ggplot2::geom_segment(ggplot2::aes(x = .data$CI.L, xend = .data$CI.R,
        yend = stats::reorder(.data$feature, .data$logFC))) + ggplot2::geom_point(size = 2) +
      ggplot2::scale_colour_manual(values = c("#74838a", "#a63350")) +
      ggplot2::labs(title = paste(name, contrast, "\u2014 effect estimates across biological replicates"),
        x = "Difference on the declared workflow scale (moderated 95% CI)", y = NULL, colour = "BH FDR < 0.05") + .wf_plot_theme(),
      caption = "Top 12 features by adjusted P value are a display subset. Inference and BH correction use the complete modality/contrast family. Confidence intervals are pointwise and do not account for selecting these features.")
    genes <- head(top$feature, 4)
    points <- do.call(rbind, lapply(genes, function(g) data.frame(feature = g,
      value = as.numeric(t$expression[g, ]), replicate = t$units$replicate, group = t$units$group)))
    p <- ggplot2::ggplot(points, ggplot2::aes(.data$group, .data$value, group = .data$replicate))
    if (isTRUE(t$contrast$paired)) p <- p + ggplot2::geom_line(colour = "#9eafb8", linewidth = .5)
    plots[[paste0("replicate_units_", name, "_", contrast)]] <- .wf_figure(p +
      ggplot2::geom_point(ggplot2::aes(colour = .data$group), size = 2) + ggplot2::facet_wrap(~feature, scales = "free_y") +
      ggplot2::labs(title = paste(name, contrast, "\u2014 inspect individual replicate means"),
        x = NULL, y = "Mean workflow value", colour = "Group") + .wf_plot_theme(),
      caption = "Each point is an equally weighted biological-replicate mean. Connecting lines indicate complete matched pairs; no lines are used for independent groups.")
  }
  for (name in names(result$analysis)) {
    p <- result$analysis[[name]]$pathways
    regional <- p$regional
    if (!is.null(regional) && nrow(regional)) {
      figure <- plotRegionalPathways(regional, num_display = 10, text_size = 10,
        database = result$resources[[name]]$pathway_index$raw_resources, verbose = FALSE)
      plots[[paste0("pathway_regional_", name)]] <- .wf_figure(figure,
        width = 14, height = 7,
        caption = paste(name, "- findRegionalPathways() uses complete identity-level effect rankings, with no foreground/FDR prefilter. Bubble size is the number of tested members. BH correction covers the complete eligible pathway family within each region. These competitive feature-set tests do not establish reproducibility across specimens."))
    }
    if (!is.null(p$geseca) && nrow(p$geseca)) {
      g <- head(as.data.frame(p$geseca)[order(p$geseca$padj), ], 12)
      g$name <- result$resources[[name]]$pathway_index$metadata$pathwayName[
        match(g$pathway, result$resources[[name]]$pathway_index$metadata$pathwayRampId)]
      if ("pctVar" %in% names(g)) plots[[paste0("pathway_geseca_", name)]] <- .wf_figure(
        ggplot2::ggplot(g, ggplot2::aes(.data$pctVar, stats::reorder(.data$name, .data$pctVar), colour = .data$padj)) +
          ggplot2::geom_point(size = 3) + ggplot2::scale_colour_viridis_c() + .wf_plot_theme() +
          ggplot2::labs(title = paste(name, "- coordinated pathway variation"), x = "GESECA explained variance", y = NULL, colour = "BH FDR"),
        width = 13, height = 6, caption = "runRAMPGeseca() uses the same measured identity assay and pathway index. This is feature-set co-regulation across observations, not a test of between-patient DE or biochemical pathway activity.")
    }
  }
  plots
}
