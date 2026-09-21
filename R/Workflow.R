# Internal helpers deliberately share the native container and pathway APIs.
.wf_named <- function(x, what) {
  if (!is.list(x) || !length(x) || is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x)))
    stop(what, " must be a non-empty, uniquely named list.", call. = FALSE)
  invisible(x)
}

.wf_number <- function(x, name, minimum = 0, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < minimum ||
      (integer && x != floor(x))) stop("Invalid ", name, ".", call. = FALSE)
  x
}

.wf_load <- function(x) {
  origin <- list(kind = "in-memory")
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop("Input RDS does not exist: ", x, call. = FALSE)
    origin <- list(kind = "RDS", path = normalizePath(x), md5 = unname(tools::md5sum(x)))
    x <- readRDS(x)
  } else if (is.list(x) && !is.null(x$resource)) {
    if (!requireNamespace("SpaMTPData", quietly = TRUE))
      stop("Install SpaMTPData to acquire registered resources.", call. = FALSE)
    if (is.null(x$version) || identical(x$version, "latest"))
      stop("Pin a SpaMTPData version in resource inputs.", call. = FALSE)
    allowed <- c("resource", "version", "local_dir", "offline")
    if (length(setdiff(names(x), allowed))) stop("Unknown resource input fields.", call. = FALSE)
    origin <- c(list(kind = "SpaMTPData"), x)
    if (is.null(x$offline)) x$offline <- TRUE
    x <- do.call(SpaMTPData::spaMTPData, x)
  }
  if (!methods::is(x, "SpatialExperiment"))
    stop("Workflow inputs must be SpatialExperiment objects or RDS/resource specifications. ",
      "Convert and bin raw spectra first with an explicit mass resolution.", call. = FALSE)
  origin$metadata <- S4Vectors::metadata(x)[intersect(
    c("SpaMTPData", "native_resource", "example", "spamtp_mapping", "spatial_alignment",
      "image_provenance", "registered_pairing", "input_exclusions"),
    names(S4Vectors::metadata(x)))]
  list(object = x, origin = origin)
}

.wf_spec <- function(x) {
  if (!is.list(x) || is.null(x$type) || is.null(x$layer))
    stop("Each modality needs explicit type and layer fields.", call. = FALSE)
  allowed <- c("type", "layer", "normalization", "scale_factor", "min_detected",
    "species", "annotation", "pathways")
  if (length(setdiff(names(x), allowed))) stop("Unknown modality specification fields: ",
    paste(setdiff(names(x), allowed), collapse = ", "), call. = FALSE)
  x$type <- match.arg(x$type, c("metabolomics", "transcriptomics", "proteomics", "other"))
  if (!is.character(x$layer) || length(x$layer) != 1L || !nzchar(x$layer))
    stop("layer must name one input assay.", call. = FALSE)
  x$normalization <- match.arg(x$normalization %||% "auto",
    c("auto", "library_log", "log1p", "none"))
  if (x$normalization == "auto") x$normalization <- if (x$layer != "counts") "none" else
    if (x$type %in% c("metabolomics", "transcriptomics")) "library_log" else
      if (x$type == "proteomics") "log1p" else "none"
  x$scale_factor <- .wf_number(x$scale_factor %||% 10000, "scale_factor", .Machine$double.eps)
  x$min_detected <- .wf_number(x$min_detected %||% 1, "min_detected", 0, TRUE)
  x
}

.wf_matrix <- function(object, spec) {
  E <- .assayData(object, layer = spec$layer)
  if (nrow(E) < 1L || ncol(E) < 1L || is.null(rownames(E)) ||
      is.null(colnames(E)) || anyNA(rownames(E)) || anyNA(colnames(E)) ||
      any(!nzchar(rownames(E))) || any(!nzchar(colnames(E))) ||
      anyDuplicated(rownames(E)) || anyDuplicated(colnames(E)))
    stop("Inputs need non-empty, unique feature and observation identifiers.", call. = FALSE)
  # Check in blocks, including delayed matrices, without densifying the full input.
  for (i in split(seq_len(nrow(E)), ceiling(seq_len(nrow(E)) / 500))) {
    values <- as.matrix(E[i, , drop = FALSE])
    if (any(!is.finite(values))) stop("Input layer contains non-finite values.", call. = FALSE)
    if (spec$normalization != "none" && any(values < 0))
      stop("Log normalization requires non-negative input; choose normalization='none' for transformed data.", call. = FALSE)
  }
  E
}

.wf_qc <- function(E, name, spec) {
  detected <- as.numeric(Matrix::colSums(E != 0))
  total <- as.numeric(Matrix::colSums(E))
  data.frame(modality = name, observation = colnames(E), total = total,
    detected = detected, pass = detected >= spec$min_detected,
    layer = spec$layer, stringsAsFactors = FALSE)
}

.wf_normalize <- function(E, spec) {
  if (spec$normalization == "none") return(E)
  if (spec$normalization == "library_log") {
    e <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = E))
    e <- normalizeSMData(e, normalisation.type = "LogNormalize",
      scale.factor = spec$scale_factor, verbose = FALSE)
    # Native normalization uses natural logs; the workflow declares log2 units.
    return(SummarizedExperiment::assay(e, "logcounts") / log(2))
  }
  log1p(E) / log(2)
}

.wf_map <- function(source, target, spec) {
  if (!is.list(spec) || is.null(spec$method))
    stop("Each moving modality needs a mapping method and geometry.", call. = FALSE)
  src <- .nativeCoordinates(source)
  dst <- .nativeCoordinates(target)
  if (spec$method == "nearest") {
    .wf_number(spec$max_distance, "mapping max_distance", .Machine$double.eps)
    if (length(setdiff(names(spec), c("method", "max_distance"))))
      stop("Unknown nearest mapping settings.", call. = FALSE)
    from <- rep(NA_integer_, nrow(dst)); distance <- rep(NA_real_, nrow(dst))
    for (sample in unique(dst$sample_id)) {
      si <- which(src$sample_id == sample); ti <- which(dst$sample_id == sample)
      if (!length(si)) next
      nn <- FNN::get.knnx(as.matrix(src[si, c("x", "y")]),
        as.matrix(dst[ti, c("x", "y")]), k = 1)
      distance[ti] <- nn$nn.dist[, 1]
      keep <- nn$nn.dist[, 1] <= spec$max_distance
      from[ti[keep]] <- si[nn$nn.index[keep, 1]]
    }
    hits <- which(!is.na(from))
    weights <- Matrix::sparseMatrix(i = from[hits], j = hits, x = 1,
      dims = c(nrow(src), nrow(dst)), dimnames = list(src$cell, dst$cell))
    matches <- lapply(from, function(i) if (is.na(i)) integer() else i)
  } else if (spec$method == "pixel") {
    allowed <- c("method", "width", "target_radius", "target_polygons", "overlap_threshold")
    if (length(setdiff(names(spec), allowed))) stop("Unknown pixel mapping settings.", call. = FALSE)
    .wf_number(spec$width, "pixel width", .Machine$double.eps)
    threshold <- .wf_number(spec$overlap_threshold %||% 0.2, "overlap_threshold")
    if (threshold > 1) stop("overlap_threshold must be <= 1.", call. = FALSE)
    mapped <- .mappingWeights(src, dst, spec$width,
      is.null(spec$target_radius) && is.null(spec$target_polygons),
      spec$target_radius, spec$target_polygons, threshold)
    weights <- mapped$weights; matches <- mapped$matches
    distance <- vapply(seq_along(matches), function(i) {
      s <- matches[[i]]
      if (!length(s)) return(NA_real_)
      min(sqrt((src$x[s] - dst$x[i])^2 + (src$y[s] - dst$y[i])^2))
    }, numeric(1))
  } else stop("Mapping method must be nearest or pixel.", call. = FALSE)
  list(weights = weights, settings = spec,
    observations = data.frame(observation = dst$cell, sample_id = dst$sample_id,
      mapped = lengths(matches) > 0, n_source = lengths(matches), distance = distance),
    matches = stats::setNames(lapply(matches, function(i) src$cell[i]), dst$cell))
}

.wf_pca <- function(E, npcs, max_features) {
  variance <- as.numeric(Matrix::rowMeans(E * E) - Matrix::rowMeans(E)^2)
  selected <- which(is.finite(variance) & variance > 1e-12)
  selected <- head(selected[order(variance[selected], decreasing = TRUE)], max_features)
  if (length(selected) < 2L || ncol(E) < 3L)
    return(list(status = "skipped", reason = "Fewer than two variable features or three observations."))
  x <- SingleCellExperiment::SingleCellExperiment(assays = list(workflow = E))
  x <- runMetabolicPCA(x, npcs = npcs, slot = "workflow",
    features = rownames(E)[selected], scale = TRUE, verbose = FALSE)
  scores <- SingleCellExperiment::reducedDim(x, "pca")
  list(status = "completed", scores = scores, features = rownames(E)[selected],
    percent_variance = attr(scores, "percentVar"), feature_scale = TRUE,
    provenance = S4Vectors::metadata(x)$reduction_inputs$pca)
}

# Compatibility adapter: the analysis itself is owned by the native DE API.
.wf_contrasts <- function(E, metadata, group, replicate, contrasts) {
  if (is.null(group)) return(list())
  x <- SingleCellExperiment::SingleCellExperiment(
    assays = list(workflow = E), colData = S4Vectors::DataFrame(metadata))
  findAllDEMs(x, ident = group, assay = "main", slot = "workflow",
    method = "replicate", replicate = replicate, contrasts = contrasts)
}

#' Run a reproducible spatial multi-omics workflow and HTML report
#'
#' Acquires native inputs, audits QC, aligns and maps independent modalities,
#' normalizes explicitly selected layers, computes per-modality PCA and an
#' equal-weight joint representation, and produces descriptive associations,
#' optional replicate-level contrasts, annotations and shared-index pathways.
#' See the End_to_End_Workflow vignette for the input and configuration schema.
#' @param input A SpatialExperiment, its RDS path, a pinned SpaMTPData resource
#'   specification, or a named list of independent inputs of these forms.
#' @param modalities Named list of specifications with type and layer, plus
#'   optional normalization, scale_factor, min_detected, species, annotation,
#'   and pathways. Paired input names are main followed by selected altExp names.
#' @param output_dir Optional new directory for report, tables and result RDS.
#' @param reference Name of the reference modality for independent inputs.
#' @param alignment Named settings for each moving modality: method='identity'
#'   explicitly declares registered coordinates; otherwise arguments for
#'   applySpatialAlignment (e.g. method='affine', alignment or landmarks).
#' @param mapping Named per-moving-modality specifications: method='nearest'
#'   with max_distance, or method='pixel' with width and optional target_radius,
#'   target_polygons and overlap_threshold, in aligned coordinate units.
#' @param group Optional colData field for descriptive summaries and plot colour.
#' @param replicate Optional colData field identifying biological replicates.
#' @param contrasts Named list of numerator, denominator and paired settings.
#'   Tests require at least three independent biological replicates per group.
#' @param observation_unit Description of what a column represents.
#' @param coordinate_units Description of reference coordinate units; required
#'   when independent inputs are mapped.
#' @param npcs,max_features Maximum PCs and variable features per modality.
#' @param clusters Optional explicit k for exploratory k-means in joint space,
#'   or in PCA space when exactly one modality is selected.
#'   NULL does not choose a biological partition automatically.
#' @param association_features Maximum regional-marker/variable-feature anchors
#'   for each modality pair. Targets are the other modality's PCA-selected
#'   features. Raw and covariate-residual correlations have no pixel P values.
#' @param seed Random seed, restored on exit.
#' @param title Report title.
#' @param regions Region metadata column, defaulting to group. NULL can use
#'   requested workflow clusters when group is absent.
#' @param native Named native spatial settings; see analyzeSpaMTPRegions.
#' @param structure Native representation/evaluation settings: spatial, primary
#'   (pca or spatial), neighbors, lambda, platform, external reference metadata,
#'   k_grid, stability_blocks and umap. Graph PCA uses all retained observations
#'   and the same selected, scaled features as ordinary PCA. Clustering stability
#'   conditions on full-data embeddings; it is not independent-specimen validation.
#' @return A spamtp_workflow list containing object, QC, alignment/mapping audit,
#'   analysis, settings, stage log and session information. output_dir also
#'   receives a self-contained report.html, CSV tables and workflow.rds.
#' @details Only observations passing QC in every selected modality and having
#'   mappings in every moving modality enter joint analysis. Excluded positions
#'   remain in the audit. Raw input assays are preserved. Automatic normalization
#'   uses log2(1 + 10000 * value / column sum) for metabolomics/transcriptomics
#'   counts, log2(1 + value) for protein counts, and preserves other layers.
#'   This is a configurable baseline, not a batch correction or platform-specific
#'   raw-signal preprocessing method. Spatial mapping can reuse source pixels;
#'   correlations and clusters are exploratory. Contrasts average workflow
#'   values within biological replicate and group, not random pixel pools.
#' @export
#' @importFrom utils capture.output head str
#' @examples
#' utils::str(formals(runSpaMTPWorkflow))
runSpaMTPWorkflow <- function(input, modalities, output_dir = NULL,
    reference = NULL, alignment = list(), mapping = list(), group = NULL,
    replicate = NULL, contrasts = list(), observation_unit = "spatial observation",
    coordinate_units = NULL, npcs = 15, max_features = 2000, clusters = NULL,
    association_features = 20, seed = 1, title = "SpaMTP multi-omics report",
    regions = group, native = list(), structure = list()) {
  .wf_named(modalities, "modalities")
  modalities <- lapply(modalities, .wf_spec)
  .wf_number(npcs, "npcs", 1, TRUE); .wf_number(max_features, "max_features", 2, TRUE)
  .wf_number(association_features, "association_features", 0, TRUE)
  .wf_number(seed, "seed", 0, TRUE)
  withr::local_seed(seed)
  if (length(contrasts) && is.null(group))
    stop("Contrasts require a group field.", call. = FALSE)
  warnings <- character()
  capture <- function(expr) withCallingHandlers(expr, warning = function(w)
    warnings <<- unique(c(warnings, conditionMessage(w))))
  if (!is.null(output_dir) && dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
    stop("output_dir must be new or empty; existing reports are preserved.", call. = FALSE)
  stages <- data.frame(stage = character(), status = character(), detail = character())
  record <- function(stage, status, detail) {
    stages[nrow(stages) + 1L, ] <<- list(stage, status, detail)
  }
  independent <- is.list(input) && is.null(input$resource)
  origins <- list(); align_audit <- list(); map_audit <- list(); sources <- list()
  if (independent) {
    .wf_named(input, "independent inputs")
    if (!setequal(names(input), names(modalities)) || length(reference) != 1L || !reference %in% names(input))
      stop("Independent inputs and modality names must match; choose a reference.", call. = FALSE)
    if (is.null(coordinate_units) || !nzchar(coordinate_units))
      stop("Declare coordinate_units for independent input mapping.", call. = FALSE)
    moving <- setdiff(names(input), reference)
    if (!setequal(names(alignment), moving) || !setequal(names(mapping), moving))
      stop("Provide alignment and mapping settings for exactly the moving modalities.", call. = FALSE)
    for (name in names(modalities)) {
      acquired <- .wf_load(input[[name]])
      sources[[name]] <- acquired$object; origins[[name]] <- acquired$origin
    }
  } else {
    if (length(alignment) || length(mapping)) stop("Alignment/mapping settings require independent inputs.", call. = FALSE)
    acquired <- .wf_load(input); parent <- acquired$object
    if (names(modalities)[1] != "main") stop("Paired modalities must start with main.", call. = FALSE)
    reference <- "main"
    origins$paired <- acquired$origin
    for (name in names(modalities)) sources[[name]] <- .experimentForAssay(parent, name)
  }
  modalities <- modalities[c(reference, setdiff(names(modalities), reference))]
  if (independent && "main" %in% setdiff(names(modalities), reference))
    stop("The name main is reserved for the reference modality.", call. = FALSE)
  record("acquisition", "completed", paste(length(sources), "modalities loaded; input metadata and resource versions retained."))
  qc <- list(); matrices <- list()
  for (name in names(modalities)) {
    matrices[[name]] <- .wf_matrix(sources[[name]], modalities[[name]])
    if (!is.null(modalities[[name]]$species))
      S4Vectors::metadata(sources[[name]])$organism <- modalities[[name]]$species
    qc[[name]] <- .wf_qc(matrices[[name]], name, modalities[[name]])
  }
  if (independent) {
    for (name in names(sources)) {
      sources[[name]] <- sources[[name]][, qc[[name]]$pass, drop = FALSE]
      matrices[[name]] <- matrices[[name]][, qc[[name]]$pass, drop = FALSE]
      if (!ncol(sources[[name]])) stop("No observations pass QC in ", name, call. = FALSE)
    }
    parent <- sources[[reference]]
    keep <- rep(TRUE, ncol(parent))
    for (name in setdiff(names(sources), reference)) {
      cfg <- alignment[[name]]
      if (is.null(cfg$method)) stop("Declare alignment settings for ", name, call. = FALSE)
      before <- .nativeCoordinates(sources[[name]])
      if (identical(cfg$method, "identity")) {
        if (length(setdiff(names(cfg), "method"))) stop("identity alignment takes no transform arguments.", call. = FALSE)
        detail <- list(method = "identity", declaration = "User declared coordinates already registered.")
      } else {
        if (any(names(cfg) %in% c("SM.data", "ST.data", "return", "store")))
          stop("Alignment container/return arguments are managed by the workflow.", call. = FALSE)
        .singleAlignmentSample(sources[[name]]); .singleAlignmentSample(parent)
        fit <- do.call(applySpatialAlignment,
          c(list(SM.data = sources[[name]], ST.data = parent, return = "result"), cfg))
        sources[[name]] <- fit$object
        detail <- fit[setdiff(names(fit), "object")]
      }
      align_audit[[name]] <- list(before = before, after = .nativeCoordinates(sources[[name]]), result = detail)
      mapped <- .wf_map(sources[[name]], parent, mapping[[name]])
      map_audit[[name]] <- mapped
      matrices[[name]] <- matrices[[name]] %*% mapped$weights
      colnames(matrices[[name]]) <- colnames(parent)
      keep <- keep & mapped$observations$mapped
    }
    record("alignment_mapping", "completed", paste(sum(keep), "of", length(keep),
      "QC-passing reference observations have all modalities; unmatched positions excluded from analysis."))
  } else {
    keep <- Reduce(`&`, lapply(qc, `[[`, "pass"))
    if (any(vapply(sources, function(x) !identical(colnames(x), colnames(parent)), logical(1))))
      stop("Paired modality observation IDs/order must match.", call. = FALSE)
    record("alignment_mapping", "provided", "Pairing supplied in altExp; no new registration validation was performed.")
  }
  retained <- data.frame(observation = colnames(parent), included = keep)
  if (sum(keep) < 3L) stop("Fewer than three observations have complete modalities after QC/mapping.", call. = FALSE)
  object <- parent[, keep, drop = FALSE]
  # Do not carry reductions/derived altExps that no longer describe this workflow.
  SingleCellExperiment::altExps(object) <- S4Vectors::SimpleList()
  SingleCellExperiment::reducedDims(object) <- S4Vectors::SimpleList()
  analysis <- list(); pcas <- list(); values <- list()
  assay_names <- stats::setNames(ifelse(names(modalities) == reference, "main", names(modalities)), names(modalities))
  for (name in names(modalities)) {
    spec <- modalities[[name]]
    E <- matrices[[name]][, keep, drop = FALSE]
    normalized <- .wf_normalize(E, spec); dimnames(normalized) <- dimnames(E)
    values[[name]] <- normalized
    if (name == reference) {
      experiment <- object
    } else if (!independent) {
      experiment <- sources[[name]][, keep, drop = FALSE]
    } else {
      experiment <- SingleCellExperiment::SingleCellExperiment(
        assays = stats::setNames(list(E), spec$layer),
        rowData = SummarizedExperiment::rowData(sources[[name]]),
        metadata = S4Vectors::metadata(sources[[name]]))
    }
    SummarizedExperiment::assay(experiment, "workflow") <- normalized
    if (!is.null(spec$species)) S4Vectors::metadata(experiment)$organism <- spec$species
    pcas[[name]] <- .wf_pca(normalized, npcs, max_features)
    if (pcas[[name]]$status == "completed") {
      SingleCellExperiment::reducedDim(experiment, paste0("workflow_", name)) <- pcas[[name]]$scores
    }
    if (name == reference) object <- experiment else
      SingleCellExperiment::altExp(object, name) <- experiment
    if (pcas[[name]]$status == "completed")
      SingleCellExperiment::reducedDim(object, paste0("workflow_", name)) <- pcas[[name]]$scores
    analysis[[name]] <- list(pca = pcas[[name]],
      comparisons = capture(.wf_contrasts(normalized, as.data.frame(SummarizedExperiment::colData(object)), group, replicate, contrasts)))
  }
  record("qc_normalization", "completed", paste(sum(keep), "observations retained across all selected modalities; raw layers preserved."))
  record("single_modality", "completed", "Variable-feature PCA and configured group summaries computed; constants excluded from PCA only.")
  joint <- NULL
  if (length(pcas) >= 2L && all(vapply(pcas, function(x) x$status == "completed", logical(1)))) {
    order <- c(reference, setdiff(names(modalities), reference))
    object <- multiOmicIntegration(object,
      modalities = unname(assay_names[order]),
      reduction.list = as.list(paste0("workflow_", order)),
      dims.list = lapply(pcas[order], function(x) seq_len(ncol(x$scores))),
      return.intermediate = TRUE)
    embedding <- SingleCellExperiment::reducedDim(object, "integrated")
    joint <- list(embedding = embedding,
      display = stats::prcomp(embedding, rank. = 2)$x[, 1:2, drop = FALSE],
      method = "Equal-weight concatenated standardized modality PCs; display is a PCA of the entire joint representation. No batch correction.")
    record("integration", "completed", joint$method)
  } else record("integration", "skipped", "Requires at least two modalities with estimable PCA.")
  clustering <- NULL
  if (!is.null(clusters)) {
    .wf_number(clusters, "clusters", 2, TRUE)
    cluster_input <- if (!is.null(joint)) joint$embedding else
      if (length(pcas) == 1L && identical(pcas[[1]]$status, "completed")) pcas[[1]]$scores else NULL
    if (is.null(cluster_input)) stop("Clustering requires an estimable joint embedding or a single-modality PCA.", call. = FALSE)
    if (clusters >= nrow(unique(cluster_input))) stop("clusters must be smaller than the number of distinct observations.", call. = FALSE)
    labels <- .wf_native_cluster(cluster_input, clusters, seed)
    object$workflow_cluster <- factor(labels)
    if (!is.null(joint)) joint$clusters <- labels
    clustering <- list(labels = labels, centers = clusters,
      input = if (!is.null(joint)) "joint embedding" else paste(names(pcas)[1], "PCA"),
      method = "k-means; requested number of exploratory clusters, not validated anatomical regions")
    record("clustering", "completed", paste(clustering$method, "using", clustering$input))
  }
  associations <- list()
  record("associations", if (length(associations)) "completed" else "skipped",
    "Descriptive Pearson correlations of selected variable features; spatial dependence and reused pixels preclude pixel-based significance claims.")
  for (name in names(modalities)) {
    spec <- modalities[[name]]; assay <- assay_names[[name]]
    if (!is.null(spec$annotation)) {
      if (spec$type != "metabolomics" || is.null(spec$annotation$index))
        stop("Workflow annotation requires a metabolomics modality and an explicit MZ index.", call. = FALSE)
      args <- spec$annotation
      if (any(names(args) %in% c("data", "assay", "return.only.annotated", "save.intermediate", "filepath")))
        stop("Annotation input/output arguments are managed by the workflow.", call. = FALSE)
      object <- capture(do.call(annotateSM, c(list(data = object, assay = assay,
        return.only.annotated = FALSE, save.intermediate = TRUE), args)))
      analysis[[name]]$annotation <- .featureMetadata(object, assay)
      record(paste0("annotation_", name), "completed", "Accurate-mass candidates with configured index; identities remain putative.")
    } else record(paste0("annotation_", name), "not_requested", "Existing feature metadata retained; no new mass annotation configured.")
    if (!is.null(spec$pathways)) {
      cfg <- spec$pathways
      if (is.null(cfg$index) || length(setdiff(names(cfg), c("index", "database", "min_size", "max_size", "foreground", "duplicate_genes", "regional", "geseca", "network", "annotation_source", "metabolite_ambiguity"))))
        stop("pathways needs index and optional min_size/max_size/foreground/duplicate_genes.", call. = FALSE)
      if (spec$type == "metabolomics") {
        identity_assay <- paste0("workflow_ids_", name)
        object <- capture(createPathwayAssay(object, analyte_type = "metabolites", assay = assay,
          slot = "workflow", new_assay = identity_assay, pathway_index = cfg$index, database = cfg$database,
          annotation_source = cfg$annotation_source %||% "current",
          metabolite_ambiguity = cfg$metabolite_ambiguity %||% "exclude", verbose = FALSE))
        target <- paste0("workflow_pathways_", name)
        context <- .pathway_score_context(object, identity_assay, "workflow", NULL, cfg$index,
          "auto", NULL, NULL, "latest", NULL, spec$species %||% "custom", "latest", "auto", NULL,
          "error", cfg$min_size %||% 5, cfg$max_size %||% 500)
        if (!any(context$coverage$eligible)) {
          analysis[[name]]$pathways <- list(status = "skipped", identity_assay = identity_assay,
            reason = "No compound pathway meets measured member limits after the declared ambiguity policy.",
            coverage = context$coverage, mapping = context$audit$inputs, provenance = context$audit$provenance,
            compound_mapping = S4Vectors::metadata(SingleCellExperiment::altExp(object, identity_assay))$pathway_mapping)
          record(paste0("pathways_", name), "skipped", "No compound pathway meets measured member limits after the declared ambiguity policy.")
          next
        }
        object <- capture(createPathwayObject(object, assay = identity_assay, slot = "workflow",
          new.assay = target, pathway_index = cfg$index, organism = spec$species %||% "custom",
          min_path_size = cfg$min_size %||% 5, max_path_size = cfg$max_size %||% 500))
        audit <- S4Vectors::metadata(SingleCellExperiment::altExp(object, target))$pathway_mapping
        analysis[[name]]$pathways <- list(assay = target, identity_assay = identity_assay,
          coverage = audit$coverage, mapping = audit$inputs, provenance = audit$provenance,
          compound_mapping = S4Vectors::metadata(SingleCellExperiment::altExp(object, identity_assay))$pathway_mapping)
        record(paste0("pathways_", name), "completed", "Explicit compound identities, configured ambiguity policy, measured members and original features retained.")
        next
      }
      if (spec$type != "transcriptomics") stop("Pathways require genes or explicitly annotated metabolites.", call. = FALSE)
      if (is.null(spec$species)) stop("Declare species before gene pathway mapping.", call. = FALSE)
      target <- paste0("workflow_pathways_", name)
      scored <- capture(tryCatch(createPathwayObject(object, assay = assay, slot = "workflow", new.assay = target,
        pathway_index = cfg$index, organism = spec$species,
        min_path_size = cfg$min_size %||% 5, max_path_size = cfg$max_size %||% 500,
        duplicate_genes = cfg$duplicate_genes %||% "error"), error = function(e) {
          if (startsWith(conditionMessage(e), "No pathways contain")) return(NULL)
          stop(e)
        }))
      if (is.null(scored)) {
        context <- capture(.pathway_score_context(object, assay, "workflow", NULL, cfg$index,
          "auto", NULL, NULL, "latest", NULL, spec$species, "latest", "auto", NULL,
          cfg$duplicate_genes %||% "error", cfg$min_size %||% 5, cfg$max_size %||% 500))
        analysis[[name]]$pathways <- list(status = "skipped", coverage = context$coverage,
          reason = "No mapped pathways satisfy the configured measured size limits.",
          mapping = context$audit$inputs, provenance = context$audit$provenance)
        record(paste0("pathways_", name), "skipped", "No mapped pathways satisfy the configured measured size limits.")
        next
      }
      object <- scored
      p <- SingleCellExperiment::altExp(object, target)
      audit <- S4Vectors::metadata(p)$pathway_mapping
      analysis[[name]]$pathways <- list(assay = target, coverage = audit$coverage,
        mapping = audit$inputs, provenance = audit$provenance)
      if (!is.null(cfg$foreground)) {
        analysis[[name]]$pathways$enrichment <- capture(fishersPathwayAnalysis(
          list(genes = cfg$foreground), universe = list(genes = rownames(values[[name]])),
          pathway_index = cfg$index, organism = spec$species,
          min_path_size = cfg$min_size %||% 5, max_path_size = cfg$max_size %||% 500,
          verbose = FALSE))
      }
      record(paste0("pathways_", name), "completed", "Shared identity index; measured-universe coverage, used members and conflicts retained. Scores are descriptive.")
    } else record(paste0("pathways_", name), "not_requested", "No species-specific pathway index configured.")
  }
  test_count <- sum(vapply(analysis, function(x) sum(vapply(x$comparisons$tests,
    function(t) identical(t$status, "completed"), logical(1))), integer(1)))
  record("replicate_contrasts", if (test_count) "completed" else if (length(contrasts)) "skipped" else "not_requested",
    paste(test_count, "replicate-level tests completed; per-test design and skipped reasons are retained."))
  settings <- list(modalities = lapply(modalities, function(x) {
    if (!is.null(x$pathways$database)) x$pathways$database <- list(
      resources = names(x$pathways$database), fingerprint = rlang::hash(x$pathways$database))
    if (!is.null(x$pathways$index)) x$pathways$index <- x$pathways$index$provenance
    if (!is.null(x$annotation$index)) {
      idx <- x$annotation$index
      x$annotation$index <- list(fingerprint = rlang::hash(idx), polarity = idx$polarity,
        maldi_matrix = idx$maldi_matrix, compounds = nrow(idx$compounds), rules = nrow(idx$rules))
    }
    x
  }), reference = reference, independent = independent, group = group, replicate = replicate,
    contrasts = contrasts, observation_unit = observation_unit, coordinate_units = coordinate_units,
    npcs = npcs, max_features = max_features, clusters = clusters,
    association_features = association_features, seed = seed, assay_names = assay_names)
  result <- structure(list(schema_version = 1L, title = title, object = object,
    qc = do.call(rbind, qc), retained = retained, alignment = align_audit, mapping = map_audit,
    analysis = analysis, joint = joint, clustering = clustering, associations = associations, settings = settings,
    stages = stages, warnings = warnings, origins = origins,
    resources = lapply(modalities, function(x) list(mass_index = x$annotation$index,
      pathway_index = x$pathways$index, pathway_database = x$pathways$database)),
    created = format(Sys.time(), tz = "UTC", usetz = TRUE),
    session = utils::sessionInfo(), package_version = as.character(utils::packageVersion("SpaMTP"))),
    class = "spamtp_workflow")
  S4Vectors::metadata(result$object)$spamtp_workflow <- list(settings = settings, stages = stages)
  result <- capture(.wf_structure(result, structure))
  result <- capture(analyzeSpaMTPRegions(result, regions = regions, native = native))
  result <- capture(.wf_associations(result))
  result <- capture(.wf_regional_pathways(result))
  result$warnings <- warnings
  if (!is.null(output_dir)) result$report <- renderSpaMTPReport(result, output_dir)
  result
}
