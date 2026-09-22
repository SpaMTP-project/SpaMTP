installedExamplePath <- function(folder, file = NULL) {
    # Prefer the source checkout when present; never audit an older installed
    # copy while testing uninstalled scripts. Built-package checks use system.file.
    path <- testthat::test_path("..", "..", "inst", folder)
    if (!dir.exists(path)) path <- system.file(folder, package = "SpaMTP")
    if (!is.null(file) && dir.exists(path)) path <- file.path(path, file)
    path
}

installedExampleEnvironment <- function(folder, file) {
    env <- new.env(parent = globalenv())
    sys.source(installedExamplePath(folder, file), envir = env)
    env
}

installedExampleViolations <- function(code) {
    tokens <- utils::getParseData(parse(text = code, keep.source = TRUE))
    bad <- tokens$token == "'@'" |
        (tokens$token == "SYMBOL_FUNCTION_CALL" &
            tokens$text %in% c("slot", "slot<-", "load_all", "FindAllMarkers")) |
        (tokens$token == "SYMBOL_PACKAGE" &
            tokens$text %in% c("Seurat", "SeuratObject", "devtools")) |
        (tokens$token == "STR_CONST" &
            grepl("^['\"](Seurat(Object)?['\"]|/(vast|stornext|home|Users)/)",
                  tokens$text))
    tokens$text[bad]
}

test_that("installed demos and workflows cannot reintroduce legacy interfaces", {
    files <- unlist(lapply(c("examples", "workflows"), function(folder) {
        list.files(installedExamplePath(folder), pattern = "\\.[Rr]$",
            full.names = TRUE, recursive = TRUE)
    }))
    expect_gte(length(files), 3L)
    for (file in files) {
        expect_length(installedExampleViolations(readLines(file, warn = FALSE)), 0L)
    }
    # Both syntactic forms count as direct slot access, even without @.
    for (bad in c("x@images", "methods::slot(x, 'image')",
                  "slot(x, 'image') <- y", "Seurat::FindAllMarkers(x)",
                  "requireNamespace('SeuratObject')",
                  "devtools::load_all('/vast/scratch/example')")) {
        expect_gt(length(installedExampleViolations(bad)), 0L)
    }
})

test_that("sourcing installed scripts only defines functions", {
    scripts <- c(examples = "pathway_network_mouse_brain_demo.R",
        workflows = "prepare_mouse_brain.R", workflows = "run_mouse_brain.R")
    for (i in seq_along(scripts)) {
        before <- loadedNamespaces()
        env <- installedExampleEnvironment(names(scripts)[i], scripts[i])
        expect_setequal(loadedNamespaces(), before)
        expect_true(all(vapply(as.list(env), is.function, logical(1))))
        expect_false(exists("args", envir = env, inherits = FALSE))
    }
})

test_that("native entry points reject invalid CLI inputs without running analysis", {
    demo <- installedExampleEnvironment("examples", "pathway_network_mouse_brain_demo.R")
    prepare <- installedExampleEnvironment("workflows", "prepare_mouse_brain.R")
    run <- installedExampleEnvironment("workflows", "run_mouse_brain.R")
    expect_error(demo$mouseBrainNetworkMain(character()), "Usage")
    expect_error(prepare$prepareMouseBrainMain(character()), "Usage")
    expect_error(run$runMouseBrainMain(c("a", "b", "c", "other")), "Usage")
    existing <- withr::local_tempdir()
    sentinel <- file.path(existing, "keep.rds")
    saveRDS("existing output", sentinel)
    expect_error(demo$mouseBrainNetworks(NULL, NULL, existing), "new output")
    expect_error(prepare$prepareMouseBrain("absent", existing), "new output")
    expect_error(run$runMouseBrain("absent", "absent", existing, "dhb"), "new output")
    expect_identical(readRDS(sentinel), "existing output")
    expect_error(demo$mouseBrainNetworks(list(), NULL, tempfile()), "SpatialExperiment")
})

test_that("image geometry is explicit, validated and retains its precision", {
    prepare <- installedExampleEnvironment("workflows", "prepare_mouse_brain.R")
    expect_identical(prepare$mouseBrainSpatialConfig(), list())
    cfg <- list(source = "Acquisition metadata", target_radius = 188.47123456789)
    expect_identical(prepare$mouseBrainSpatialConfig(cfg), cfg)
    for (value in list(0, -1, NA_real_, Inf, "100", c(1, 2))) {
        invalid <- cfg
        invalid$target_radius <- value
        expect_error(prepare$mouseBrainSpatialConfig(invalid), "target_radius")
    }
    expect_error(prepare$mouseBrainSpatialConfig(list(target_radius = 1)), "provenance")
    expect_error(prepare$mouseBrainSpatialConfig(c(cfg, list(typo = 1))), "accepts")
    cfg$images <- list(mouse_brain_visium = list(path = "missing.png", scale_factor = 1))
    expect_error(prepare$mouseBrainSpatialConfig(cfg), "existing path")
    file <- tempfile(fileext = ".png")
    on.exit(unlink(file), add = TRUE)
    file.create(file)
    cfg$images$mouse_brain_visium$path <- file
    expect_identical(prepare$mouseBrainSpatialConfig(cfg)$images$mouse_brain_visium$path,
        normalizePath(file))
    cfg$images$mouse_brain_visium$scale_factor <- 0
    expect_error(prepare$mouseBrainSpatialConfig(cfg), "scale_factor")
    run <- installedExampleEnvironment("workflows", "run_mouse_brain.R")
    folder <- withr::local_tempdir()
    geometry <- list(pixel_width = 500, target_radius = NULL,
        coordinate_units = "original image pixels")
    saveRDS(geometry, file.path(folder, "pairing.rds"))
    expect_error(run$mouseBrainPairing(folder), "target_radius")
    geometry$target_radius <- 188.47123456789
    saveRDS(geometry, file.path(folder, "pairing.rds"))
    expect_identical(run$mouseBrainPairing(folder), geometry)
    geometry$coordinate_units <- "micrometres"
    saveRDS(geometry, file.path(folder, "pairing.rds"))
    expect_error(run$mouseBrainPairing(folder), "original image pixels")
})

test_that("native preparation preserves paired data and embeds optional optical images", {
    skip_if_not_installed("SpaMTPData", "0.99.5")
    skip_if_not(requireNamespace("SpaMTP", quietly = TRUE), "Native dependencies unavailable")
    prepare <- installedExampleEnvironment("workflows", "prepare_mouse_brain.R")
    object <- nativeFixture()
    parent <- withr::local_tempdir()
    image <- file.path(parent, "optical.png")
    png::writePNG(array(0.5, c(8, 8, 3)), image)
    seen <- list()
    testthat::local_mocked_bindings(spaMTPData = function(resource, version,
            local_dir, offline) {
        seen[[resource]] <<- list(version = version, offline = offline)
        object
    }, .package = "SpaMTPData")
    output <- file.path(parent, "prepared")
    config <- list(source = "Synthetic geometry for testing only",
        target_radius = 0.123456789012345,
        images = list(mouse_brain_visium = list(path = image, scale_factor = 0.5)))
    prepare$prepareMouseBrain(parent, output, config)
    expect_length(seen, 3L)
    expect_true(all(vapply(seen, function(x) identical(x,
        list(version = "1.1.0", offline = TRUE)), logical(1))))
    # NULL uses the public reader's configured local/cache/download resolution.
    # An explicit offline flag still forbids downloading with no directory.
    prepare$prepareMouseBrain(outputDir = file.path(parent, "online"))
    expect_true(all(vapply(seen, function(x) identical(x,
        list(version = "1.1.0", offline = FALSE)), logical(1))))
    prepare$prepareMouseBrain(outputDir = file.path(parent, "cached"), offline = TRUE)
    expect_true(all(vapply(seen, function(x) identical(x,
        list(version = "1.1.0", offline = TRUE)), logical(1))))
    visium <- readRDS(file.path(output, "visium.rds"))
    expect_identical(SummarizedExperiment::assays(visium), SummarizedExperiment::assays(object))
    expect_identical(SingleCellExperiment::altExp(visium, "transcriptome"),
        SingleCellExperiment::altExp(object, "transcriptome"))
    expect_equal(SpatialExperiment::spatialCoords(visium), SpatialExperiment::spatialCoords(object))
    expect_true(methods::validObject(visium))
    expect_equal(SpatialExperiment::imgData(visium)$sample_id, unique(visium$sample_id))
    expect_equal(SpatialExperiment::imgData(visium)$scaleFactor, 0.5)
    expect_equal(S4Vectors::metadata(visium)$image_provenance$image_md5,
        unname(tools::md5sum(image)))
    geometry <- jsonlite::read_json(file.path(output, "pairing.json"), simplifyVector = TRUE)
    expect_equal(geometry$target_radius, config$target_radius, tolerance = 1e-14)
    # A loaded image remains usable independently of its original file path.
    unlink(image)
    expect_true(length(SpatialExperiment::imgRaster(visium)) > 0L)
    expect_equal(nrow(SpatialExperiment::imgData(object)), 0L)
})

test_that("workflow readers preserve configured resources and explicit offline mode", {
    skip_if_not_installed("SpaMTPData", "0.99.5")
    run <- installedExampleEnvironment("workflows", "run_mouse_brain.R")
    demo <- installedExampleEnvironment("examples", "pathway_network_mouse_brain_demo.R")
    parent <- withr::local_tempdir()
    saveRDS(nativeFixture(), file.path(parent, "dhb.rds"))
    withr::local_options(list(SpaMTPdb.resource_dir = "configured-resources"))
    calls <- list()
    testthat::local_mocked_bindings(loadSpaMTPDatabase = function(resources,
            version, local_dir, offline) {
        calls$db <<- list(version = version, local_dir = local_dir,
            offline = offline, configured = getOption("SpaMTPdb.resource_dir"))
        stop("reader captured")
    }, .package = "SpaMTP")
    testthat::local_mocked_bindings(spaMTPData = function(resource, version,
            local_dir, offline) {
        calls$data <<- list(version = version, local_dir = local_dir, offline = offline)
        nativeFixture()
    }, .package = "SpaMTPData")
    expect_error(run$runMouseBrain(parent, outputDir = tempfile(), case = "dhb"),
        "reader captured")
    expect_identical(calls$db, list(version = "3.0.7", local_dir = NULL,
        offline = FALSE, configured = "configured-resources"))
    expect_error(run$runMouseBrain(parent, "explicit-resources", tempfile(), "dhb"),
        "reader captured")
    expect_identical(calls$db, list(version = "3.0.7", local_dir = "explicit-resources",
        offline = TRUE, configured = "explicit-resources"))
    expect_identical(getOption("SpaMTPdb.resource_dir"), "configured-resources")
    expect_error(demo$runMouseBrainNetworkDemo(outputDir = tempfile()), "reader captured")
    expect_identical(calls$data, list(version = "1.1.0", local_dir = NULL, offline = FALSE))
    expect_false(calls$db$offline)
    expect_identical(calls$db$configured, "configured-resources")
    expect_error(demo$runMouseBrainNetworkDemo(outputDir = tempfile(), offline = TRUE),
        "reader captured")
    expect_true(calls$data$offline)
    expect_true(calls$db$offline)
    expect_error(demo$runMouseBrainNetworkDemo("native", "database", tempfile()),
        "reader captured")
    expect_identical(calls$data$local_dir, "native")
    expect_identical(calls$db$local_dir, "database")
    expect_true(calls$data$offline)
    expect_true(calls$db$offline)
})

test_that("the installed network demo uses native effects and preserves paired RNA", {
    skip_if_not(requireNamespace("SpaMTP", quietly = TRUE), "Native dependencies unavailable")
    demo <- installedExampleEnvironment("examples", "pathway_network_mouse_brain_demo.R")
    object <- nativeFixture()
    ids <- paste0("RAMP_C_", 1:3)
    db <- list(source_df = data.frame(rampId = ids, commonName = letters[1:3],
        sourceId = paste0("hmdb:", 1:3)),
        chem_props = data.frame(ramp_id = ids, common_name = letters[1:3],
            chem_source_id = paste0("hmdb:", 1:3)),
        pathway = data.frame(pathwayRampId = "P1", pathwayName = "Fixture pathway",
            sourceId = "WP1", type = "WikiPathways"),
        analytehaspathway = data.frame(rampId = ids[1:2], pathwayRampId = "P1"),
        ramp_db_metadata = list(ramp_version = "fixture"),
        ramp_wikipathway = list(list(id = "WP1", title = "Fixture pathway",
            mixedEdges = data.frame(src = ids[1:2], dest = ids[2:3],
                directed = 1L, reaction_type = 1L))))
    # Isolate annotation chemistry (tested by the annotation-engine suite), but
    # run real normalization, compound aggregation, scran, enrichment and HTML.
    annotationArgs <- NULL
    local_mocked_bindings(annotateSM = function(data, ...) {
        annotationArgs <<- list(...)
        annotations <- data.frame(mz_name = rownames(data),
            observed_mz = c(100, 200, 300), Adduct = "M+H",
            Ramp_IDs = ids, Score = .9, MassScore = .99, ChemicalScore = 1,
            IsotopeScore = NA_real_, AdductNetworkScore = NA_real_)
        S4Vectors::metadata(data)$mz_annotation <- list(results = annotations,
            metadata = list(engine = "indexed-chemical-v2", ramp_version = "fixture"))
        data
    }, .package = "SpaMTP")
    parent <- withr::local_tempdir()
    result <- demo$mouseBrainNetworks(object, db, file.path(parent, "networks"),
        ident = "region", minPathSize = 2)
    expect_false(annotationArgs$return.only.annotated)
    expect_identical(annotationArgs$maldi_matrix, "DHB")
    expect_identical(result$status, "completed")
    expect_identical(SummarizedExperiment::assay(result$object, "counts"),
        SummarizedExperiment::assay(object, "counts"))
    expect_identical(SingleCellExperiment::altExp(result$object, "transcriptome"),
        SingleCellExperiment::altExp(object, "transcriptome"))
    expect_false(any(c("p_val_adj", "P.Value", "FDR") %in% names(result$markers$DEMs)))
    expect_setequal(names(result$networks), c("leading_edge", "annotated"))
    expect_true(all(file.exists(unlist(result$networks))))
    expect_length(unique(dirname(unlist(result$networks))), 2L)
    expect_true(file.exists(file.path(parent, "networks", "network_analysis.rds")))
    empty <- demo$mouseBrainNetworks(object, db, file.path(parent, "no-pathways"),
        ident = "region", minPathSize = 4)
    expect_identical(empty$status, "no_eligible_pathways")
    expect_length(empty$networks, 0L)
    expect_true(file.exists(file.path(parent, "no-pathways", "network_analysis.rds")))
})
