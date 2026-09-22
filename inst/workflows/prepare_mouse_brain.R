# Prepare native mouse-brain inputs, optionally attaching independent images.
# Rscript prepare_mouse_brain.R NATIVE_1.1.0_DIR NEW_OUTPUT_DIR [SPATIAL_CONFIG.json]
# Sourcing this file only defines functions. No historical object is read.

mouseBrainSpatialConfig <- function(config = NULL) {
    if (is.null(config)) return(list())
    if (is.character(config) && length(config) == 1L) {
        if (!file.exists(config)) stop("Spatial configuration file not found.")
        config <- jsonlite::read_json(config, simplifyVector = FALSE)
    }
    allowed <- c("source", "target_radius", "images")
    if (!is.list(config) || is.data.frame(config) ||
        is.null(names(config)) || anyDuplicated(names(config)) ||
        any(!names(config) %in% allowed)) {
        stop("Spatial configuration accepts source, target_radius and images only.")
    }
    if (!is.character(config$source) || length(config$source) != 1L ||
        is.na(config$source) || !nzchar(trimws(config$source))) {
        stop("Provide source provenance for the image scale and/or spot radius.")
    }
    positive <- function(x) is.numeric(x) && length(x) == 1L &&
        is.finite(x) && x > 0
    if (!is.null(config$target_radius) && !positive(config$target_radius)) {
        stop("target_radius must be positive, in original Visium image pixels.")
    }
    if (!is.null(config$images)) {
        resources <- c("mouse_brain_visium", "mouse_brain_dhb_striatum")
        if (!is.list(config$images) || is.null(names(config$images)) ||
            anyDuplicated(names(config$images)) ||
            any(!names(config$images) %in% resources)) {
            stop("images must be named by the Visium or DHB resource identifier.")
        }
        for (name in names(config$images)) {
            image <- config$images[[name]]
            if (!is.list(image) || !setequal(names(image), c("path", "scale_factor")) ||
                anyDuplicated(names(image)) ||
                !is.character(image$path) || length(image$path) != 1L ||
                is.na(image$path) || !file.exists(image$path) ||
                dir.exists(image$path) || !positive(image$scale_factor)) {
                stop("Each image needs an existing path and a positive scale_factor.")
            }
            config$images[[name]]$path <- normalizePath(image$path, mustWork = TRUE)
        }
    }
    config
}

# NULL directories use public readers (verified local files, cache or download).
# Supplying a directory retains the offline command-line recipe's default.
prepareMouseBrain <- function(nativeDir = NULL, outputDir, spatialConfig = NULL,
                              offline = !is.null(nativeDir)) {
    if (length(outputDir) != 1L || is.na(outputDir) || !nzchar(outputDir) ||
        file.exists(outputDir)) {
        stop("Supply a new output directory; existing paths are not overwritten.")
    }
    config <- mouseBrainSpatialConfig(spatialConfig)
    for (package in c("SpaMTP", "SpaMTPData")) {
        if (!requireNamespace(package, quietly = TRUE)) {
            stop("Install the native package and its dependencies: ", package)
        }
    }
    if (utils::packageVersion("SpaMTP") < "0.99.9" ||
        utils::packageVersion("SpaMTPData") < "0.99.5") {
        stop("Use SpaMTP >= 0.99.9 and SpaMTPData >= 0.99.5 for this recipe.")
    }
    resources <- c(visium = "mouse_brain_visium", fmp10 = "mouse_brain_fmp10",
        dhb = "mouse_brain_dhb_striatum")
    objects <- lapply(resources, function(resource) {
        x <- SpaMTPData::spaMTPData(resource, version = "1.1.0",
            local_dir = nativeDir, offline = offline)
        if (!methods::is(x, "SpatialExperiment") ||
            length(unique(x$sample_id)) != 1L) {
            stop("Expected a single-sample native resource: ", resource)
        }
        x
    })
    # This pairing applies only to these registered, transformed source inputs.
    # Resource names are not a pairing rule for arbitrary mouse-brain datasets.
    specimen <- "mouse_brain_FMP10_Visium"
    for (name in c("visium", "fmp10")) {
        x <- objects[[name]]
        x$sample_id <- rep(specimen, ncol(x))
        objects[[name]] <- x
    }
    xy <- SpatialExperiment::spatialCoords(objects$fmp10)
    if (!all(c("x", "y") %in% colnames(xy)) || any(!is.finite(xy))) {
        stop("FMP10 needs finite native x/y coordinates.")
    }
    xy <- unique(xy[, c("x", "y"), drop = FALSE])
    if (nrow(xy) < 2L) stop("FMP10 needs at least two distinct pixel centres.")
    spacing <- FNN::get.knn(xy, k = 1)$nn.dist[, 1]
    width <- stats::median(spacing[is.finite(spacing) & spacing > 0])
    if (!is.finite(width) || width <= 0) stop("Cannot determine the MSI pixel width.")
    pairing <- list(specimen = specimen,
        evidence = "Registered transformed FMP10/Visium pair from the published mouse-brain example",
        alignment = "already registered in the pinned source resources",
        source_versions = list(native = "1.1.0"),
        coordinate_units = "original image pixels", pixel_width = width,
        target_radius = config$target_radius,
        width_method = "Median nearest-neighbour distance between distinct native MSI pixel centres",
        target_radius_method = if (is.null(config$target_radius))
            "Not supplied; required only for FMP10/Visium pixel-overlap mapping" else
            "Explicit radius in original Visium image pixels from spatial configuration",
        spatial_source = config$source,
        resources = lapply(objects, function(x) S4Vectors::metadata(x)$SpaMTPData))
    for (name in c("visium", "fmp10")) {
        S4Vectors::metadata(objects[[name]])$registered_pairing <- pairing
    }
    for (name in names(objects)) {
        resource <- resources[[name]]
        image <- config$images[[resource]]
        if (is.null(image)) next
        objects[[name]] <- SpaMTP::addSpatialImage(objects[[name]],
            imageSource = image$path, scaleFactor = image$scale_factor,
            sampleId = unique(objects[[name]]$sample_id), imageId = "histology",
            load = TRUE)
        S4Vectors::metadata(objects[[name]])$image_provenance <- list(
            source = config$source, resource = resource,
            file_name = basename(image$path),
            image_md5 = unname(tools::md5sum(image$path)),
            scale_factor = image$scale_factor,
            coordinate_units = "original image pixels")
    }
    for (x in objects) methods::validObject(x)
    if (!dir.create(outputDir, recursive = TRUE)) stop("Cannot create output directory.")
    for (name in names(objects)) {
        saveRDS(objects[[name]], file.path(outputDir, paste0(name, ".rds")))
    }
    saveRDS(pairing, file.path(outputDir, "pairing.rds"))
    jsonlite::write_json(pairing, file.path(outputDir, "pairing.json"),
        pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null")
    if (is.null(config$target_radius)) {
        message("DHB is ready. FMP10/Visium mapping additionally requires an explicit target_radius.")
    }
    invisible(pairing)
}

prepareMouseBrainMain <- function(args = commandArgs(trailingOnly = TRUE)) {
    if (!length(args) %in% c(2L, 3L)) {
        stop(paste("Usage: Rscript prepare_mouse_brain.R",
            "NATIVE_1.1.0_DIR NEW_OUTPUT_DIR [SPATIAL_CONFIG.json]"))
    }
    config <- if (length(args) == 3L) args[3] else NULL
    prepareMouseBrain(args[1], args[2], config)
}

if (sys.nframe() == 0L) prepareMouseBrainMain()
