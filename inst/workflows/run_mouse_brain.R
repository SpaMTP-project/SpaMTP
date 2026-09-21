# Rscript run_mouse_brain.R PREPARED_INPUT_DIR DB_3.0.7_DIR NEW_OUTPUT_DIR CASE
# CASE is fmp10 or dhb. Run prepare_mouse_brain.R first.
# Sourcing this file only defines functions.

mouseBrainPairing <- function(inputDir) {
    path <- file.path(inputDir, "pairing.rds")
    if (!file.exists(path)) stop("Run prepare_mouse_brain.R first.")
    geometry <- readRDS(path)
    if (!is.list(geometry) || is.data.frame(geometry)) {
        stop("Pairing geometry must be a named list from native preparation.")
    }
    for (field in c("pixel_width", "target_radius")) {
        value <- geometry[[field]]
        if (!is.numeric(value) || length(value) != 1L ||
            !is.finite(value) || value <= 0) {
            stop("FMP10 mapping requires a positive ", field,
                ". Supply target_radius (original image pixels) in the ",
                "preparation spatial configuration; it is not inferred.")
        }
    }
    if (!identical(geometry$coordinate_units, "original image pixels")) {
        stop("Pairing geometry must declare original image pixels; re-run preparation.")
    }
    geometry
}

runMouseBrain <- function(inputDir, databaseDir, outputDir,
                          case = c("fmp10", "dhb")) {
    case <- match.arg(case)
    if (length(outputDir) != 1L || is.na(outputDir) || !nzchar(outputDir) ||
        file.exists(outputDir)) stop("Supply a new output directory.")
    geometry <- if (case == "fmp10") mouseBrainPairing(inputDir) else NULL
    inputs <- if (case == "fmp10") c("visium.rds", "fmp10.rds") else "dhb.rds"
    if (!all(file.exists(file.path(inputDir, inputs)))) {
        stop("Prepared native inputs are missing; run prepare_mouse_brain.R first.")
    }
    if (!requireNamespace("SpaMTP", quietly = TRUE)) {
        stop("Install SpaMTP and its native dependencies.")
    }
    if (utils::packageVersion("SpaMTP") < "0.99.9") {
        stop("Use SpaMTP >= 0.99.9 for this recipe.")
    }
    withr::local_options(list(SpaMTPdb.resource_dir = databaseDir))
    db <- SpaMTP::loadSpaMTPDatabase(c("chem_props", "source_df",
        "analytehaspathway", "pathway", "ramp_db_metadata",
        "ramp_wikipathway", "ramp_reactome", "ramp_kegg", "ramp_hmdb"),
        version = "3.0.7", local_dir = databaseDir, offline = TRUE)
    index <- SpaMTP::buildPathwayIndex(db,
        gene_mapping = "ramp", organism = "Mus musculus")
    matrixName <- if (case == "fmp10") "FMP-10" else "DHB"
    massIndex <- SpaMTP::buildMZAnnotationIndex(db$chem_props,
        polarity = "positive", maldi_matrix = matrixName)
    ms <- list(type = "metabolomics", layer = "counts", species = "Mus musculus",
        annotation = list(index = massIndex, ppm_error = 5, min_score = 0,
            verbose = FALSE),
        pathways = list(index = index, database = db, min_size = 3,
            max_size = 100, annotation_source = "current",
            metabolite_ambiguity = "exclude", geseca = TRUE, network = TRUE))
    rna <- list(type = "transcriptomics", layer = "counts",
        species = "Mus musculus")
    # The HGNC index is human-specific. Mouse gene identities remain unchanged.
    if (case == "dhb") {
        result <- SpaMTP::runSpaMTPWorkflow(file.path(inputDir, "dhb.rds"),
            list(main = ms, SPT = rna), group = "lesion", regions = "RegionLoupe",
            npcs = 15, max_features = 2000, clusters = 6, association_features = 12,
            coordinate_units = "original image pixels",
            observation_unit = "matched spatial observation from one DHB/Visium specimen",
            structure = list(primary = "spatial", reference = "RegionLoupe",
                k_grid = c(4, 6, 8, 10)),
            title = "Mouse brain DHB: native spatial multi-omics analysis")
    } else {
        reference <- readRDS(file.path(inputDir, "visium.rds"))
        if (!methods::is(reference, "SpatialExperiment") ||
            !"RegionLoupe" %in% colnames(SummarizedExperiment::colData(reference))) {
            stop("visium.rds must be a SpatialExperiment with RegionLoupe labels.")
        }
        labelled <- !is.na(reference$RegionLoupe) &
            nzchar(trimws(as.character(reference$RegionLoupe)))
        S4Vectors::metadata(reference)$input_exclusions <- data.frame(
            observation = colnames(reference)[!labelled],
            reason = "Outside the annotated analysis area")
        reference <- reference[, labelled]
        result <- SpaMTP::runSpaMTPWorkflow(
            list(RNA = reference, MS = file.path(inputDir, "fmp10.rds")),
            list(RNA = rna, MS = ms), reference = "RNA",
            group = "lesion", regions = "RegionLoupe",
            npcs = 15, max_features = 2000, clusters = 10, association_features = 12,
            coordinate_units = geometry$coordinate_units,
            observation_unit = "Visium spot with overlapping registered MSI pixels",
            alignment = list(MS = list(method = "identity")),
            mapping = list(MS = list(method = "pixel", width = geometry$pixel_width,
                target_radius = geometry$target_radius, overlap_threshold = 0.3)),
            structure = list(primary = "spatial", reference = "RegionLoupe",
                k_grid = c(6, 8, 10, 12)),
            title = "Mouse brain FMP10 / Visium: native registered analysis")
    }
    if (!dir.create(outputDir, recursive = TRUE)) stop("Cannot create output directory.")
    result$report <- SpaMTP::renderSpaMTPReport(result, file.path(outputDir, "report"))
    saveRDS(result, file.path(outputDir, "analysis.rds"))
    print(result$structure$metrics)
    message("Report: ", result$report)
    invisible(result)
}

runMouseBrainMain <- function(args = commandArgs(trailingOnly = TRUE)) {
    if (length(args) != 4L || !args[4] %in% c("fmp10", "dhb")) {
        stop(paste("Usage: Rscript run_mouse_brain.R",
            "PREPARED_INPUT_DIR DB_3.0.7_DIR NEW_OUTPUT_DIR fmp10|dhb"))
    }
    runMouseBrain(args[1], args[2], args[3], args[4])
}

if (sys.nframe() == 0L) runMouseBrainMain()
