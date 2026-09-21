# Native, offline Mouse Brain/DHB pathway-network example.
# Rscript pathway_network_mouse_brain_demo.R NATIVE_1.1.0_DIR DB_3.0.7_DIR NEW_OUTPUT_DIR [REGION_COLUMN]
# Alternatively, source this file and call runMouseBrainNetworkDemo(). Sourcing
# only defines functions: it does not download data, attach packages or analyse.

mouseBrainNetworks <- function(object, database, outputDir,
                               ident = "RegionLoupe", minPathSize = 3L,
                               maxPathSize = 100L) {
    if (length(outputDir) != 1L || is.na(outputDir) || !nzchar(outputDir) ||
        file.exists(outputDir)) {
        stop("Supply a new output directory; existing paths are not overwritten.")
    }
    if (!methods::is(object, "SpatialExperiment")) {
        stop("Supply a native SpatialExperiment, not a historical object.")
    }
    if (length(unique(object$sample_id)) != 1L) {
        stop("This single-specimen example requires exactly one sample.")
    }
    if (length(ident) != 1L || is.na(ident) ||
        !ident %in% colnames(SummarizedExperiment::colData(object))) {
        stop("ident must name a colData field containing the region labels.")
    }
    labels <- as.character(SummarizedExperiment::colData(object)[[ident]])
    keep <- !is.na(labels) & nzchar(trimws(labels))
    excluded <- colnames(object)[!keep]
    object <- object[, keep]
    if (length(unique(labels[keep])) < 2L) {
        stop("At least two labelled regions are required.")
    }
    counts <- SummarizedExperiment::assay(object, "counts")
    mz <- SummarizedExperiment::rowData(object)$mz
    if (!is.numeric(mz) || length(mz) != nrow(object) ||
        any(!is.finite(mz)) || any(mz <= 0)) {
        stop("rowData(object)$mz must contain positive finite masses.")
    }
    spectrum <- data.frame(mz = mz, intensity = Matrix::rowMeans(counts))
    # Recompute scored annotations; do not reuse the archive's mass-only cache.
    object <- SpaMTP::annotateSM(object, db = database$chem_props,
        assay = "main", raw.mz.column = "mz", ppm_error = 5,
        polarity = "positive", maldi_matrix = "DHB",
        ms1_spectrum = spectrum, min_score = 0,
        return.only.annotated = FALSE, save.intermediate = TRUE)
    object <- SpaMTP::normalizeSMData(object,
        normalisation.type = "LogNormalize", verbose = FALSE)
    SummarizedExperiment::assay(object, "workflow") <-
        SummarizedExperiment::assay(object, "logcounts") / log(2)

    # Mouse genes are retained in altExp(), but are not mapped to human HGNC.
    # This demonstration tests compound sets, not mouse gene pathways.
    index <- SpaMTP::buildPathwayIndex(database,
        gene_mapping = "ramp", organism = "Mus musculus")
    object <- SpaMTP::createPathwayAssay(object, assay = "main",
        slot = "workflow", new_assay = "compound", database = database,
        pathway_index = index, organism = "Mus musculus",
        annotation_source = "current", annotation_score_threshold = 0.05,
        metabolite_ambiguity = "exclude", verbose = FALSE)
    compounds <- SingleCellExperiment::altExp(object, "compound")
    result <- list(object = object, excluded_observations = excluded,
        pathway_index = index$provenance,
        compound_mapping = S4Vectors::metadata(compounds)$pathway_mapping,
        annotation = SpaMTP::annotationInfo(object), networks = list(),
        settings = list(region = ident, analyte_types = "metabolites",
            organism = "Mus musculus", min_path_size = minPathSize,
            max_path_size = maxPathSize, annotation_score_threshold = 0.05,
            metabolite_ambiguity = "exclude", ppm = 5, maldi_matrix = "DHB",
            interpretation = paste("Single-specimen descriptive region effects;",
                "competitive pathway P values are not biological-replicate tests.",
                "RNA is retained but is not mapped to a human gene reference.")))
    if (nrow(compounds) >= 2L) {
        markers <- SpaMTP::findAllDEMs(object, ident = ident,
            assay = "compound", slot = "workflow", method = "markers",
            spatial_blocks = 0, verbose = FALSE)
        result$markers <- markers
        if (identical(markers$status, "completed")) {
            # Retain every measured identity in each rank vector. Region
            # characterisation does not manufacture pixel-level DE P values.
            ranks <- vapply(colnames(markers$expression), function(region) {
                tab <- markers$DEMs[markers$DEMs$cluster == region, ]
                tab$logFC[match(rownames(markers$expression), tab$gene)]
            }, numeric(nrow(markers$expression)))
            dimnames(ranks) <- dimnames(markers$expression)
            result$regional <- withr::with_seed(1234,
                SpaMTP::findRegionalPathways(object, ident = ident,
                    ranks = list(metabolites = ranks), pathway_index = index,
                    min_path_size = minPathSize, max_path_size = maxPathSize,
                    organism = "Mus musculus", verbose = FALSE))
        }
    }
    if (!dir.create(outputDir, recursive = TRUE)) {
        stop("Cannot create output directory: ", outputDir)
    }
    eligible <- !is.null(result$regional) && nrow(result$regional) > 0L &&
        any(is.finite(result$regional$NES))
    result$status <- if (eligible) "completed" else "no_eligible_pathways"
    if (eligible) {
        for (mode in c("leading_edge", "annotated")) {
            folder <- file.path(outputDir, mode)
            dir.create(folder)
            result$networks[[mode]] <- SpaMTP::pathwayNetworkPlots(object,
                ident = ident, regpathway = result$regional,
                DE.list = list(metabolites = result$markers$DEMs),
                SM_assay = "compound", SM_slot = "workflow",
                analyte_types = "metabolites", metabolite_detection = mode,
                database = database, pathway_index = index,
                organism = "Mus musculus", image = NULL, path = folder,
                top_n_pathways = 4, max_nodes = 350,
                max_spatial_points = 10000, verbose = FALSE)
        }
    } else {
        message("No eligible compound pathways; annotation and coverage audits are saved.")
    }
    result$session <- utils::sessionInfo()
    saveRDS(result, file.path(outputDir, "network_analysis.rds"))
    invisible(result)
}

runMouseBrainNetworkDemo <- function(nativeDir, databaseDir, outputDir,
                                    ident = "RegionLoupe") {
    if (file.exists(outputDir)) stop("Supply a new output directory.")
    for (package in c("SpaMTP", "SpaMTPData")) {
        if (!requireNamespace(package, quietly = TRUE)) {
            stop("Install the native package and its dependencies: ", package)
        }
    }
    if (utils::packageVersion("SpaMTP") < "0.99.9" ||
        utils::packageVersion("SpaMTPData") < "0.99.5") {
        stop("Use SpaMTP >= 0.99.9 and SpaMTPData >= 0.99.5 for this example.")
    }
    withr::local_options(list(SpaMTPdb.resource_dir = databaseDir))
    object <- SpaMTPData::spaMTPData("mouse_brain_dhb_striatum",
        version = "1.1.0", local_dir = nativeDir, offline = TRUE)
    database <- SpaMTP::loadSpaMTPDatabase(c("chem_props", "source_df",
        "analytehaspathway", "pathway", "ramp_db_metadata",
        "ramp_wikipathway", "ramp_reactome", "ramp_kegg", "ramp_hmdb"),
        version = "3.0.7", local_dir = databaseDir, offline = TRUE)
    mouseBrainNetworks(object, database, outputDir, ident = ident)
}

mouseBrainNetworkMain <- function(args = commandArgs(trailingOnly = TRUE)) {
    if (!length(args) %in% c(3L, 4L)) {
        stop(paste("Usage: Rscript pathway_network_mouse_brain_demo.R",
            "NATIVE_1.1.0_DIR DB_3.0.7_DIR NEW_OUTPUT_DIR [REGION_COLUMN]"))
    }
    ident <- if (length(args) == 4L) args[4] else "RegionLoupe"
    runMouseBrainNetworkDemo(args[1], args[2], args[3], ident)
}

if (sys.nframe() == 0L) mouseBrainNetworkMain()
