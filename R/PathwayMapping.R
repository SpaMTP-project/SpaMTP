.pathway_links <- function(table) {
  required <- c("rampId", "pathwayRampId")
  if (is.null(table)) return(data.frame(rampId = character(), pathwayRampId = character()))
  if (!all(required %in% names(table))) stop("Pathway links require rampId and pathwayRampId.", call. = FALSE)
  table <- as.data.frame(table)[, required, drop = FALSE]
  table[] <- lapply(table, as.character)
  table <- unique(table[!is.na(table$rampId) & nzchar(table$rampId) &
                         !is.na(table$pathwayRampId) & nzchar(table$pathwayRampId), ])
  table <- table[order(table$pathwayRampId, table$rampId, method = "radix"), ]
  rownames(table) <- NULL
  table
}

.pathway_resource_fingerprint <- function(resources) {
  meta <- resources$pathway
  fields <- intersect(c("pathwayRampId", "pathwayName", "sourceId", "type", "pathwayCategory"), names(meta))
  meta <- as.data.frame(meta)[, fields, drop = FALSE]
  meta[] <- lapply(meta, as.character)
  if (nrow(meta)) meta <- unique(meta[do.call(order, c(unname(meta), list(na.last = TRUE, method = "radix"))), , drop = FALSE])
  rownames(meta) <- NULL
  rlang::hash(list(links = .pathway_links(resources$analytehaspathway), metadata = meta))
}

.pathway_validate_index <- function(index, resources = NULL, gene_index = NULL) {
  if (!inherits(index, "spamtp_pathway_index") || !identical(index$provenance$schema_version, 1L)) {
    stop("pathway_index must be built by buildPathwayIndex().", call. = FALSE)
  }
  if (!is.null(resources)) {
    if (!identical(.pathway_resource_fingerprint(resources), index$provenance$membership_fingerprint)) {
      stop("pathway_index does not match pathway memberships or metadata; rebuild it.", call. = FALSE)
    }
    if (!is.null(index$gene_index) && !is.null(resources$source_df)) {
      .gene_check_source_index(resources$source_df, index$gene_index)
    }
  }
  if (!is.null(gene_index) && !identical(gene_index, index$gene_index)) {
    stop("gene_index and pathway_index use different gene identities.", call. = FALSE)
  }
  invisible(TRUE)
}

.pathway_database <- function(database, pathway_index) {
  if (is.null(pathway_index)) return(database)
  .pathway_validate_index(pathway_index)
  if (!is.null(database)) return(database)
  database <- pathway_index$raw_resources
  if (is.null(database$source_df)) database$source_df <- data.frame(
    rampId = character(), sourceId = character())
  database
}

.pathway_context <- function(resources, database = resources,
    gene_mapping = "auto", gene_reference = NULL, gene_index = NULL,
    gene_reference_version = "latest", gene_reference_local_dir = NULL,
    organism = "Homo sapiens", pathway_index = NULL) {
  gene_mapping <- match.arg(gene_mapping, c("auto", "hgnc", "ramp"))
  if (!is.null(pathway_index)) {
    .pathway_validate_index(pathway_index, resources, gene_index)
    if (!is.null(pathway_index$gene_index) &&
        !identical(organism, pathway_index$gene_index$provenance$organism)) {
      stop("pathway_index organism does not match organism.", call. = FALSE)
    }
    if (identical(gene_mapping, "ramp") && !is.null(pathway_index$gene_index)) {
      stop("gene_mapping = 'ramp' conflicts with the HGNC pathway_index.", call. = FALSE)
    }
    if (identical(gene_mapping, "hgnc") && is.null(pathway_index$gene_index)) {
      stop("gene_mapping = 'hgnc' requires an HGNC pathway_index.", call. = FALSE)
    }
    resources[names(pathway_index$resources)] <- pathway_index$resources
    return(list(resources = resources,
                index = pathway_index$gene_index, pathway_index = pathway_index))
  }
  raw <- resources
  if (is.null(resources$source_df)) {
    if (!is.null(gene_index) || !is.null(gene_reference) || identical(gene_mapping, "hgnc")) {
      stop("Gene harmonization requires source_df or a prebuilt pathway_index.", call. = FALSE)
    }
    view <- list(resources = resources, index = NULL)
  } else {
    view <- .gene_pathway_view(resources, database, gene_mapping, gene_reference,
      gene_index, gene_reference_version, gene_reference_local_dir, organism)
  }
  links <- .pathway_links(view$resources$analytehaspathway)
  raw_links <- .pathway_links(raw$analytehaspathway)
  ids <- sort(unique(c(as.character(raw$pathway$pathwayRampId), raw_links$pathwayRampId)))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  sets <- function(table) {
    answer <- split(table$rampId, table$pathwayRampId)
    answer <- answer[match(ids, names(answer))]
    answer[lengths(answer) == 0L] <- list(character())
    stats::setNames(answer, ids)
  }
  meta <- as.data.frame(raw$pathway)
  if (!"pathwayRampId" %in% names(meta)) meta <- data.frame(pathwayRampId = character())
  meta <- meta[match(ids, meta$pathwayRampId), , drop = FALSE]
  meta$pathwayRampId <- ids
  if (!"pathwayName" %in% names(meta)) meta$pathwayName <- ids
  absent <- is.na(meta$pathwayName) | !nzchar(meta$pathwayName)
  meta$pathwayName[absent] <- ids[absent]
  rownames(meta) <- NULL
  conflicts <- if (is.null(view$index)) character() else {
    view$index$nodes$rampId[view$index$nodes$status == "conflicting_record"]
  }
  index <- structure(list(
    members = sets(links), raw_members = sets(raw_links), metadata = meta,
    conflicting_members = sets(raw_links[raw_links$rampId %in% conflicts, ]),
    gene_index = view$index,
    resources = view$resources[intersect(c("source_df", "analyte", "analytehaspathway", "pathway"), names(view$resources))],
    raw_resources = raw[intersect(c("source_df", "analyte", "analytehaspathway", "pathway"), names(raw))],
    provenance = list(schema_version = 1L,
      membership_fingerprint = .pathway_resource_fingerprint(raw),
      gene_mapping = if (is.null(view$index)) "ramp" else "hgnc",
      gene = view$index$provenance,
      resources = lapply(raw[intersect(c("source_df", "analytehaspathway", "pathway"), names(raw))],
                         attr, which = "spamtp_database"))
  ), class = "spamtp_pathway_index")
  list(resources = view$resources, index = view$index, pathway_index = index)
}

#' Build a shared, auditable pathway membership index
#'
#' Unites pathway memberships of non-conflicting RaMP records belonging to one
#' HGNC gene. Retains raw memberships and excluded conflicting records. Reuse
#' this index in enrichment, expression aggregation, scoring and named pathway
#' plots to obtain the same identities and membership sets. Pathways are keyed
#' by pathwayRampId; identical display names never merge distinct pathways.
#'
#' @param database A resource bundle with analytehaspathway and pathway, and
#'   source_df for gene harmonization. NULL retrieves these from SpaMTPdb.
#' @inheritParams fishersPathwayAnalysis
#' @inheritParams buildGeneMappingIndex
#' @return A spamtp_pathway_index containing canonical and raw member sets,
#'   conflicting members, metadata, gene_index and resource fingerprints.
#'   Database size means unique analyte identities after reconciliation and
#'   conflict exclusion, before intersecting with measured features.
#' @export
#' @examples
#' db <- list(analytehaspathway = data.frame(rampId = "RAMP_C_1",
#'            pathwayRampId = "P1"), pathway = data.frame(pathwayRampId = "P1"))
#' buildPathwayIndex(db)$members
buildPathwayIndex <- function(database = NULL, gene_mapping = c("auto", "hgnc", "ramp"),
    gene_reference = NULL, gene_index = NULL, gene_reference_version = "latest",
    gene_reference_local_dir = NULL, organism = "Homo sapiens",
    database_version = "latest", database_source = c("auto", "spamtpdb"),
    database_local_dir = NULL) {
  needed <- c("analytehaspathway", "pathway")
  if (is.null(database) || "source_df" %in% names(database)) needed <- c(needed, "source_df")
  resources <- .spamtp_db_bundle(needed, database, version = database_version,
    source = match.arg(database_source), local_dir = database_local_dir)
  if (!is.null(database$analyte)) resources$analyte <- database$analyte
  if (is.null(resources$analyte)) resources$analyte <- data.frame(
    rampId = unique(c(resources$source_df$rampId, resources$analytehaspathway$rampId)))
  .pathway_context(resources, database, match.arg(gene_mapping), gene_reference,
    gene_index, gene_reference_version, gene_reference_local_dir, organism)$pathway_index
}

.pathway_coverage <- function(index, measured_ids, used_ids = measured_ids,
    foreground_ids = character(), types = c("G", "C"), min_size = 1, max_size = Inf) {
  select <- function(ids) ids[grepl(paste0("^RAMP_(", paste(types, collapse = "|"), ")_"), ids)]
  # Custom feature-set fixtures may use identifiers other than RaMP.
  if (is.null(types)) select <- function(ids) ids
  members <- lapply(index$members, select)
  raw <- lapply(index$raw_members, select)
  # Index members are already unique character vectors. Subsetting avoids
  # repeated S4 set-operation dispatch for large registries with empty sets.
  measured <- lapply(members, function(ids) ids[ids %in% measured_ids])
  used <- lapply(measured, function(ids) ids[ids %in% used_ids])
  foreground <- lapply(measured, function(ids) ids[ids %in% foreground_ids])
  coverage <- index$metadata
  coverage$database_size <- lengths(members)
  coverage$database_raw_size <- lengths(raw)
  coverage$measured_size <- lengths(measured)
  coverage$coverage_fraction <- ifelse(lengths(members) > 0, lengths(measured) / lengths(members), NA_real_)
  coverage$used_size <- lengths(used)
  coverage$eligible <- lengths(used) >= min_size & lengths(used) <= max_size & lengths(used) > 0
  coverage$database_members <- I(unname(members))
  coverage$measured_members <- I(unname(measured))
  coverage$used_members <- I(unname(used))
  coverage$foreground_members <- I(unname(foreground))
  coverage$excluded_conflict_records <- I(unname(lapply(index$conflicting_members, select)))
  coverage$excluded_conflict_count <- lengths(coverage$excluded_conflict_records)
  attr(coverage, "pathway_index") <- index$provenance
  coverage
}

.pathway_sets_from_coverage <- function(coverage, eligible_only = TRUE) {
  rows <- if (eligible_only) which(coverage$eligible) else seq_len(nrow(coverage))
  stats::setNames(unclass(coverage$used_members[rows]), coverage$pathwayRampId[rows])
}

.pathway_select <- function(index, pathways) {
  if (!is.character(pathways) || !length(pathways) || anyNA(pathways)) {
    stop("pathways must contain pathway IDs or unambiguous names.", call. = FALSE)
  }
  vapply(pathways, function(query) {
    if (query %in% index$metadata$pathwayRampId) return(query)
    hit <- index$metadata$pathwayRampId[tolower(index$metadata$pathwayName) == tolower(query)]
    if (length(hit) != 1L) stop("Unknown or ambiguous pathway name: ", query,
      ". Supply an exact pathwayRampId.", call. = FALSE)
    hit
  }, character(1), USE.NAMES = FALSE)
}

.pathway_expression <- function(E, index, duplicate_genes = c("error", "mean", "sum"),
    feature_ids = rownames(E), verbose = TRUE) {
  duplicate_genes <- match.arg(duplicate_genes)
  if (is.null(feature_ids) || length(feature_ids) != nrow(E) || anyNA(feature_ids) ||
      any(!nzchar(feature_ids))) stop("Expression features need non-missing identifiers.", call. = FALSE)
  ids <- as.character(feature_ids)
  canonical <- ids
  genes <- !grepl("^RAMP_C_", ids, ignore.case = TRUE)
  audit <- NULL
  if (!is.null(index$gene_index) && any(genes)) {
    audit <- .map_gene_index(ids[genes], index$gene_index, verbose = verbose)
    canonical[genes] <- audit$canonical_ramp_id
  }
  if (is.null(index$gene_index) && any(genes) && !is.null(index$raw_resources$source_df)) {
    source <- index$raw_resources$source_df
    if (!"sourceId" %in% names(source)) source$sourceId <- NA_character_
    mapped <- .fisher_map_ids(ids, "G", source, unique(c(source$rampId, unlist(index$members))))
    by_input <- split(mapped$rampId, mapped$input)
    canonical[genes] <- NA_character_
    for (i in which(genes)) {
      match_ids <- unique(by_input[[tolower(ids[i])]])
      if (length(match_ids) > 1L) stop("Ambiguous RaMP-only expression mapping: ", ids[i], call. = FALSE)
      if (length(match_ids)) canonical[i] <- match_ids
    }
    audit <- data.frame(input = ids[genes], canonical_ramp_id = canonical[genes],
      status = ifelse(is.na(canonical[genes]), "unmapped", "mapped_ramp_only"),
      method = "ramp_source")
    if (verbose && anyNA(canonical[genes])) warning(
      sum(is.na(canonical[genes])), " gene input(s) unresolved in RaMP-only mapping; inspect the mapping audit.",
      call. = FALSE)
  }
  keep <- !is.na(canonical)
  values <- E[keep, , drop = FALSE]
  canonical <- canonical[keep]
  if (!length(canonical)) stop("No expression features have usable identities.", call. = FALSE)
  if (any(!is.finite(values))) stop("Pathway expression requires finite values.", call. = FALSE)
  groups <- split(seq_along(canonical), factor(canonical, levels = unique(canonical)))
  choose <- lapply(groups, function(rows) {
    if (length(rows) > 1L && duplicate_genes == "error") {
      if (!all(vapply(rows[-1L], function(i) isTRUE(all.equal(
          as.numeric(values[i, ]), as.numeric(values[rows[1L], ]), tolerance = 0)), logical(1)))) {
        stop("Multiple expression features map to ", canonical[rows[1L]],
             "; choose duplicate_genes = 'mean' or 'sum' after reviewing the assay scale.", call. = FALSE)
      }
      rows <- rows[1L]
    }
    rows
  })
  j <- unlist(choose, use.names = FALSE)
  i <- rep(seq_along(choose), lengths(choose))
  weight <- if (duplicate_genes == "mean") rep(1 / lengths(choose), lengths(choose)) else rep(1, length(j))
  W <- Matrix::sparseMatrix(i = i, j = j, x = weight, dims = c(length(groups), nrow(values)))
  result <- W %*% values
  dimnames(result) <- list(names(groups), colnames(E))
  list(expression = result, mapping = audit,
       original_features = stats::setNames(lapply(groups, function(rows) ids[keep][rows]), names(groups)),
       duplicate_genes = duplicate_genes)
}

.pathway_gene_ids <- function(ids, index) {
  if (!is.null(index$gene_index)) {
    mapped <- .map_gene_index(ids, index$gene_index, verbose = FALSE)$canonical_ramp_id
    return(unique(mapped[!is.na(mapped)]))
  }
  source <- index$raw_resources$source_df
  if (is.null(source)) return(unique(ids))
  if (!"sourceId" %in% names(source)) source$sourceId <- NA_character_
  unique(.fisher_map_ids(ids, "G", source, unique(c(source$rampId, unlist(index$members))))$rampId)
}

.pathway_score_matrix <- function(E, sets, standardize = TRUE) {
  scores <- lapply(sets, function(ids) {
    ids <- intersect(ids, rownames(E))
    if (!length(ids)) return(rep(NA_real_, ncol(E)))
    value <- as.numeric(Matrix::colSums(E[ids, , drop = FALSE])) / sqrt(length(ids))
    value <- value - mean(value)
    if (standardize) {
      deviation <- stats::sd(value)
      value <- if (is.finite(deviation) && deviation > 0) value / deviation else rep(0, length(value))
    }
    value
  })
  if (!length(scores)) return(matrix(numeric(), 0L, ncol(E), dimnames = list(character(), colnames(E))))
  answer <- do.call(rbind, scores)
  dimnames(answer) <- list(names(sets), colnames(E))
  answer
}

.pathway_score_context <- function(object, assay, slot, database, pathway_index,
    gene_mapping, gene_reference, gene_index, gene_reference_version,
    gene_reference_local_dir, organism, database_version, database_source,
    database_local_dir, duplicate_genes, min_path_size = 1, max_path_size = Inf) {
  # Validate the assay before touching external resources.
  E <- .assayData(object, assay, slot)
  database <- .pathway_database(database, pathway_index)
  if (is.null(pathway_index)) {
    pathway_index <- buildPathwayIndex(database, gene_mapping, gene_reference,
      gene_index, gene_reference_version, gene_reference_local_dir, organism,
      database_version, database_source, database_local_dir)
  } else {
    .pathway_context(database, database, gene_mapping, gene_reference, gene_index,
      gene_reference_version, gene_reference_local_dir, organism, pathway_index)
  }
  if (!is.null(pathway_index$gene_index)) .gene_check_experiment_species(object, assay, organism)
  prior <- S4Vectors::metadata(.experimentForAssay(object, assay))$pathway_mapping$provenance
  if (!is.null(prior) && identical(prior$gene_mapping, "hgnc") &&
      !identical(prior, pathway_index$provenance)) {
    stop("The input assay was mapped with a different pathway index. Rebuild it from original features.", call. = FALSE)
  }
  mapped <- .pathway_expression(E, pathway_index, duplicate_genes)
  types <- if (all(grepl("^RAMP_[GC]_", rownames(mapped$expression)))) {
    unique(sub("^RAMP_([GC])_.*", "\\1", rownames(mapped$expression)))
  } else NULL
  if (length(min_path_size) != 1L || !is.finite(min_path_size) || min_path_size < 1 ||
      min_path_size != floor(min_path_size) || length(max_path_size) != 1L ||
      is.na(max_path_size) || max_path_size < min_path_size) {
    stop("Invalid pathway size limits.", call. = FALSE)
  }
  coverage <- .pathway_coverage(pathway_index, rownames(mapped$expression), types = types,
    min_size = min_path_size, max_size = max_path_size)
  list(index = pathway_index, expression = mapped$expression, coverage = coverage,
       audit = list(provenance = pathway_index$provenance, inputs = mapped$mapping,
         original_features = mapped$original_features, duplicate_genes = mapped$duplicate_genes,
         coverage = coverage))
}
