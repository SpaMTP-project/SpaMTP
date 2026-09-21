.wf_native_widgets <- function(result, directory, max_points) {
  html <- list(associations = character(), pathways = character())
  status <- data.frame(function_name = character(), context = character(), status = character(), detail = character())
  record <- function(fn, context, state, detail) {
    status[nrow(status) + 1L, ] <<- list(fn, context, state, detail)
  }
  xy <- SpatialExperiment::spatialCoords(result$object)
  if (!all(c("x", "y") %in% colnames(xy))) return(list(html = html, status = status))
  dir.create(file.path(directory, "widgets"), showWarnings = FALSE)
  count <- 0L
  embed <- function(path, title) {
    encoded <- jsonlite::base64_enc(readBin(path, "raw", n = file.info(path)$size))
    paste0('<details class="native-widget"><summary>', .wf_escape(title),
      '</summary><p><a href="widgets/', basename(path), '">Open standalone interactive view</a></p>',
      '<iframe loading="lazy" title="', .wf_escape(title), '" src="data:text/html;base64,',
      encoded, '" style="width:100%;height:680px;border:0"></iframe></details>')
  }
  if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) {
    for (name in names(result$associations)) {
      t <- result$associations[[name]]
      t <- t[is.finite(t$correlation), , drop = FALSE]
      if (!nrow(t)) next
      best <- t[1, ]
      for (sample in unique(as.character(result$object$sample_id))) {
        columns <- which(result$object$sample_id == sample)
        columns <- columns[unique(as.integer(round(seq(1, length(columns), length.out = min(length(columns), max_points)))))]
        object <- result$object[, columns]
        widget <- plot3DFeature(object, features = c(best$feature1, best$feature2),
          assays = unname(result$settings$assay_names[c(best$modality1, best$modality2)]),
          slots = "workflow", names = paste(c(best$modality1, best$modality2), c(best$feature1, best$feature2)),
          between.layer.height = diff(range(SpatialExperiment::spatialCoords(object)[, "x"])) * .15,
          sampleId = sample, size = 2, plot.width = 1100, plot.height = 600)
        count <- count + 1L; path <- file.path(directory, "widgets", sprintf("colocalization_%03d.html", count))
        htmlwidgets::saveWidget(widget, path, selfcontained = TRUE)
        html$associations <- c(html$associations, embed(path, paste("Rotate paired spatial layers:", name, sample)))
        record("plot3DFeature", paste(name, sample), "rendered", paste(length(columns), "displayed observations; selected by residual correlation, not independent validation"))
      }
    }
  } else record("plot3DFeature", "HTML export", "skipped", "Self-contained native Plotly export requires rmarkdown and Pandoc. Static paired analyses remain available.")
  for (name in names(result$analysis)) {
    cfg <- result$settings$modalities[[name]]$pathways
    p <- result$analysis[[name]]$pathways
    if (!isTRUE(cfg$network) || is.null(p$regional) || !nrow(p$regional)) next
    if (length(unique(result$object$sample_id)) != 1L) {
      record("pathwayNetworkPlots", name, "skipped", "The native network viewer currently requires a single specimen; no samples were silently pooled.")
      next
    }
    type <- if (result$settings$modalities[[name]]$type == "metabolomics") "metabolites" else "genes"
    folder <- file.path(directory, "widgets", paste0("network_", match(name, names(result$analysis))))
    dir.create(folder, showWarnings = FALSE)
    path <- tryCatch(pathwayNetworkPlots(result$object, ident = result$region_analysis$field,
      regpathway = p$regional, DE.list = stats::setNames(list(p$identity_markers$DEMs), type),
      path = folder, SM_assay = p$identity_assay, ST_assay = p$identity_assay,
      SM_slot = "workflow", ST_slot = "workflow", analyte_types = type,
      pathway_index = result$resources[[name]]$pathway_index,
      organism = result$settings$modalities[[name]]$species %||% "custom",
      database = result$resources[[name]]$pathway_database,
      top_n_pathways = 3, max_spatial_points = max_points, verbose = FALSE),
      error = function(e) {
        if (grepl("stored topological structure|valid network edges", conditionMessage(e))) {
          record("pathwayNetworkPlots", name, "skipped", conditionMessage(e)); return(NULL)
        }
        stop(e)
      })
    if (!is.null(path)) {
      target <- file.path(directory, "widgets", paste0("network_", match(name, names(result$analysis)), ".html"))
      file.copy(as.character(path), target)
      html$pathways <- c(html$pathways, embed(target, paste("Explore pathway members and tissue distributions:", name)))
      record("pathwayNetworkPlots", name, "rendered", "Uses the same identity assay, regional effects and membership index as enrichment. Marker effects have no pixel-based FDR.")
    }
  }
  list(html = html, status = status)
}
