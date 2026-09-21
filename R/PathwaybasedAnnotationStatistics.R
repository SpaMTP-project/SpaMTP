.annotationStatisticsInput <- function(data, mz.assay, mz.slot, pathway.assay,
                                        pathway.slot, database, threshold,
                                        corrWeight, nWeight) {
  if (!methods::is(data, "SingleCellExperiment"))
    stop("Use SingleCellExperiment/SpatialExperiment; convert Seurat input explicitly.",
         call. = FALSE)
  if (length(threshold) != 1L || !is.finite(threshold) || abs(threshold) > 1)
    stop("corr_theshold must be between -1 and 1.", call. = FALSE)
  if (length(corrWeight) != 1L || length(nWeight) != 1L ||
      any(!is.finite(c(corrWeight, nWeight))) || min(corrWeight, nWeight) < 0 ||
      corrWeight + nWeight == 0)
    stop("Weights must be finite, non-negative and not both zero.", call. = FALSE)
  expression <- .nativeExpression(data, mz.assay, mz.slot)
  pathways <- .nativeExpression(data, pathway.assay, pathway.slot)
  if (ncol(expression) < 3L)
    stop("At least three paired pixels are required.", call. = FALSE)
  if (!all(c("rampId", "pathwayRampId") %in% names(database$analytehaspathway)) ||
      !all(c("rampId", "commonName") %in% names(database$source_df)))
    stop("Database requires pathway membership and RaMP/common-name columns.", call. = FALSE)
  members <- unique(as.data.frame(database$analytehaspathway)[,
    c("rampId", "pathwayRampId")])
  members <- members[!is.na(members$rampId) & !is.na(members$pathwayRampId), ]
  sets <- lapply(split(as.character(members$pathwayRampId), members$rampId), unique)
  candidates <- .annotationStatisticsCandidates(data, mz.assay, database$source_df)
  list(expression = expression, pathways = pathways, sets = sets,
       candidates = candidates, source = database$source_df)
}

.annotationStatisticsCandidates <- function(object, assay, source) {
  object <- .pathwayAnnotationObject(object, assay)
  features <- .featureMetadata(object)
  current <- .storedData(object, "mz_annotation")
  compatible <- .storedData(object, "db_3")
  if (.annotation_has_current_schema(current$results) ||
      .annotation_has_current_schema(compatible)) {
    value <- .resolve_pathway_metabolite_annotations(object,
      annotation_source = "current", score_threshold = 0)
    # Current annotation stores may predate arbitrary feature IDs.
    if (!all(value$mz_name %in% rownames(features))) {
      masses <- .massValues(features, rownames(features))
      missing <- !value$mz_name %in% rownames(features)
      value$mz_name[missing] <- rownames(features)[
        match(value$observed_mz[missing], masses)]
    }
    value <- unique(value[, c("mz_name", "ramp_id")])
    return(split(value$ramp_id, value$mz_name))
  }
  column <- intersect(c("all_Ramp_IDs", "Ramp_IDs", "all_Isomers_IDs"), names(features))
  if (!length(column))
    stop("No RaMP candidates in current annotations or rowData.", call. = FALSE)
  column <- column[1L]
  values <- lapply(as.character(features[[column]]), function(value) {
    ids <- trimws(unlist(strsplit(value, ";", fixed = TRUE)))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (column == "all_Isomers_IDs") {
      if (!"sourceId" %in% names(source))
        stop("Legacy annotation mapping requires source_df$sourceId.", call. = FALSE)
      ids <- source$rampId[source$sourceId %in% ids]
    }
    sort(unique(as.character(ids[!is.na(ids)])))
  })
  stats::setNames(values, rownames(features))
}

.annotationStatisticsName <- function(id, source) {
  matched <- as.character(source$commonName[source$rampId %in% id])
  matched <- matched[!is.na(matched) & nzchar(matched)]
  if (!length(matched)) return(NA_character_)
  names(sort(table(matched), decreasing = TRUE))[1L]
}

.annotationZ <- function(values) {
  deviation <- stats::sd(values)
  if (!is.finite(deviation) || deviation == 0) return(rep(0, length(values)))
  (values - mean(values)) / deviation
}

.scoreAnnotationFeature <- function(feature, input, threshold, corrWeight, nWeight) {
  ids <- sort(unique(input$candidates[[feature]]))
  if (length(ids) < 2L) return(NULL)
  selected <- sort(intersect(unique(unlist(input$sets[ids])),
                              rownames(input$pathways)))
  if (!length(selected)) return(NULL)
  target <- as.numeric(input$expression[feature, ])
  if (any(!is.finite(target))) stop("m/z expression must be finite.", call. = FALSE)
  # One target and only its candidate pathways are materialized, not a
  # feature-by-pathway all-pairs matrix or a synthetic MSI object.
  correlations <- vapply(selected, function(pathway) {
    values <- as.numeric(input$pathways[pathway, ])
    if (any(!is.finite(values)))
      stop("Pathway expression must be finite.", call. = FALSE)
    if (stats::sd(target) == 0 || stats::sd(values) == 0) return(NA_real_)
    stats::cor(target, values, method = "pearson")
  }, numeric(1))
  maximum <- vapply(ids, function(id) {
    values <- correlations[intersect(input$sets[[id]], selected)]
    values <- values[is.finite(values)]
    if (!length(values)) return(NA_real_)
    values[which.max(abs(values))]
  }, numeric(1))
  number <- vapply(ids, function(id) {
    values <- correlations[intersect(input$sets[[id]], selected)]
    as.integer(sum(values > threshold, na.rm = TRUE))
  }, integer(1))
  valid <- is.finite(maximum)
  scores <- pvalue <- adjusted <- rep(NA_real_, length(ids))
  if (sum(valid) > 1L) {
    scores[valid] <- corrWeight * .annotationZ(abs(maximum[valid])) +
      nWeight * .annotationZ(number[valid])
    pvalue[valid] <- stats::pnorm(scores[valid], lower.tail = FALSE)
    adjusted[valid] <- stats::p.adjust(pvalue[valid], method = "BH")
  }
  result <- tibble::tibble(
    metabolite = unname(vapply(ids, .annotationStatisticsName, character(1),
                        source = input$source)),
    ramp_id = ids, n_sig_path = unname(number), max_cor = unname(maximum),
    z_score = scores, pval = pvalue, pval_adj = adjusted)
  result[order(-result$z_score, result$ramp_id, na.last = TRUE), ]
}

#' Rank candidate metabolite annotations using paired pathway expression
#'
#' Pearson correlations are calculated directly between a selected m/z assay
#' and pathway expression in altExp, aligned by pixel IDs. This replaces the
#' former temporary Seurat assay and synthetic-mass Cardinal conversion.
#'
#' @details For each candidate, max_cor is the signed correlation of the
#' pathway with largest absolute correlation; ties use pathway ID order.
#' n_sig_path counts finite correlations strictly greater than corr_theshold
#' (zero therefore counts positive correlations, not every pathway).
#' The score is corr_weight * scale(abs(max_cor)) +
#' n_weight * scale(n_sig_path), across candidates for one feature.
#' A constant component contributes zero. Constant expression and unmeasured
#' pathways have undefined correlations. Scores require two evaluable
#' candidates; otherwise scores and probabilities are NA.
#'
#' pval is the standard-normal upper tail of this weighted score; pval_adj
#' applies BH within each feature. These are legacy heuristic ranking measures,
#' not calibrated identification p-values: the weighted components are not
#' independent, the combined score need not be standard normal, and pixels are
#' not biological replicates. Correlation is across all supplied pixels;
#' subset samples/conditions beforehand when appropriate.
#'
#' Current scored RaMP annotations take precedence over rowData RaMP IDs.
#' Older all_Isomers_IDs columns can be mapped via source_df$sourceId.
#' Duplicate candidate IDs are counted once.
#'
#' @param mz One numeric m/z query or exact character feature ID.
#' @param data A SingleCellExperiment or SpatialExperiment with paired altExp.
#' @param mz.assay Primary experiment ("main") or alternative MSI experiment.
#' @param pathway.assay Alternative experiment containing pathway-level scores.
#' @param mz.slot Expression assay in the MSI experiment.
#' @param pathway.slot Expression assay in the pathway experiment.
#' @param corr_theshold Minimum signed pathway correlation (legacy spelling).
#' @param corr_weight,n_weight Non-negative weights, not both zero.
#' @param database Optional named list from loadSpaMTPDatabase.
#' @param database_version,database_source,database_local_dir Database selection;
#'   see loadSpaMTPDatabase.
#' @return A tibble ordered by decreasing score, containing metabolite,
#'   ramp_id, n_sig_path, max_cor, z_score, pval and pval_adj.
#' @export
#' @examples
#' utils::str(formals(calculateSingleAnnotationStatistics))
calculateSingleAnnotationStatistics <- function(
    mz, data, mz.assay = "main", pathway.assay = "pathway",
    mz.slot = "counts", pathway.slot = "pathwayScores", corr_theshold = 0,
    corr_weight = 1, n_weight = 1, database = NULL,
    database_version = "latest", database_source = c("auto", "spamtpdb"),
    database_local_dir = NULL
) {
  resources <- .spamtp_db_bundle(c("source_df", "analytehaspathway"),
    database = database, version = database_version,
    source = match.arg(database_source), local_dir = database_local_dir)
  input <- .annotationStatisticsInput(data, mz.assay, mz.slot, pathway.assay,
    pathway.slot, resources, corr_theshold, corr_weight, n_weight)
  if (length(mz) != 1L || is.na(mz)) stop("Supply one m/z query.", call. = FALSE)
  if (is.numeric(mz)) {
    if (!is.finite(mz)) stop("m/z must be finite.", call. = FALSE)
    masses <- .massValues(.featureMetadata(data, mz.assay), rownames(input$expression))
    mz <- rownames(input$expression)[which.min(abs(masses - mz))]
  }
  if (!mz %in% rownames(input$expression)) stop("Unknown m/z feature.", call. = FALSE)
  result <- .scoreAnnotationFeature(mz, input, corr_theshold, corr_weight, n_weight)
  if (is.null(result))
    stop("The feature needs at least two candidate IDs and measured pathways.", call. = FALSE)
  result
}

#' Rank annotations for all MSI features
#'
#' When pathway.assay contains RaMP analytes rather than pathway scores,
#' createPathwayObject is reused to calculate scores in a temporary altExp.
#' Set pathway.scores=TRUE to use an existing pathway-level assay directly.
#' Input containers are never modified.
#' @details Uses the score definition and missing-data rules documented in
#' [calculateSingleAnnotationStatistics()]. The reported pval and pval_adj
#' columns are heuristic ranking measures, not calibrated identification
#' p-values or biological-replicate inference.
#' @inheritParams calculateSingleAnnotationStatistics
#' @param pathway.assay Alternative experiment of RaMP analytes, or pathway
#'   scores when pathway.scores=TRUE.
#' @param pathway.slot Input expression assay in pathway.assay.
#' @param return.top Return one best evaluable candidate per feature; FALSE
#'   returns a named list of all candidate tables. Features without a score
#'   remain as NA rows in top results and NULL entries when no table is possible.
#' @param pathway.scores Whether pathway.assay already contains pathway scores.
#' @return A data.frame preserving feature order and rowData, with ranked
#'   annotation columns; or a named list when return.top=FALSE.
#' @inheritParams createPathwayAssay
#' @export
#' @examples
#' utils::str(formals(calculateAnnotationStatistics))
calculateAnnotationStatistics <- function(
    data, mz.assay = "main", pathway.assay = "pathway",
    mz.slot = "counts", pathway.slot = "counts", return.top = TRUE,
    corr_theshold = 0, corr_weight = 1, n_weight = 1,
    database = NULL, database_version = "latest",
    database_source = c("auto", "spamtpdb"), database_local_dir = NULL,
    pathway.scores = FALSE, pathway_index = NULL,
    gene_mapping = c("auto", "hgnc", "ramp"), gene_reference = NULL, gene_index = NULL,
    gene_reference_version = "latest", gene_reference_local_dir = NULL,
    organism = "Homo sapiens", duplicate_genes = c("error", "mean", "sum")
) {
  source <- match.arg(database_source)
  database <- .pathway_database(database, pathway_index)
  required <- c("source_df", "analytehaspathway", if (!pathway.scores) "pathway")
  resources <- .spamtp_db_bundle(required, database = database,
    version = database_version, source = source, local_dir = database_local_dir)
  if (!methods::is(data, "SingleCellExperiment"))
    stop("Use SingleCellExperiment/SpatialExperiment; convert Seurat input explicitly.",
         call. = FALSE)
  if (!pathway.scores) {
    .nativeExpression(data, pathway.assay, pathway.slot)
    newName <- tail(make.unique(c(SingleCellExperiment::altExpNames(data),
                                  ".annotationPathways")), 1L)
    data <- createPathwayObject(data, assay = pathway.assay, slot = pathway.slot,
      new.assay = newName, database = resources, database_version = database_version,
      database_source = source, database_local_dir = database_local_dir,
      pathway_index = pathway_index, gene_mapping = match.arg(gene_mapping),
      gene_reference = gene_reference %||% database$gene_reference, gene_index = gene_index,
      gene_reference_version = gene_reference_version, gene_reference_local_dir = gene_reference_local_dir,
      organism = organism, duplicate_genes = match.arg(duplicate_genes))
    pathway.assay <- newName
    pathway.slot <- "pathwayScores"
  }
  input <- .annotationStatisticsInput(data, mz.assay, mz.slot, pathway.assay,
    pathway.slot, resources, corr_theshold, corr_weight, n_weight)
  ids <- rownames(input$expression)
  result <- stats::setNames(lapply(ids, .scoreAnnotationFeature, input = input,
    threshold = corr_theshold, corrWeight = corr_weight, nWeight = n_weight), ids)
  audit <- S4Vectors::metadata(.experimentForAssay(data, pathway.assay))$pathway_mapping
  if (!is.null(pathway_index) && !is.null(audit) && !identical(audit$provenance, pathway_index$provenance)) {
    stop("Stored pathway scores use a different pathway_index.", call. = FALSE)
  }
  if (!return.top) {
    attr(result, "pathway_coverage") <- audit$coverage
    return(result)
  }
  empty <- data.frame(metabolite = NA_character_, ramp_id = NA_character_,
    n_sig_path = NA_integer_, max_cor = NA_real_, z_score = NA_real_,
    pval = NA_real_, pval_adj = NA_real_)
  top <- do.call(rbind, lapply(result, function(table) {
    if (is.null(table) || !any(is.finite(table$z_score))) return(empty)
    as.data.frame(table[which(is.finite(table$z_score))[1L], ])
  }))
  metadata <- .featureMetadata(data, mz.assay)
  metadata$mz_names <- ids
  metadata$mz <- .massValues(metadata, ids)
  # A repeated call must not create duplicate output column names.
  for (column in names(top)) metadata[[column]] <- top[[column]]
  rownames(metadata) <- ids
  attr(metadata, "pathway_coverage") <- audit$coverage
  metadata
}
