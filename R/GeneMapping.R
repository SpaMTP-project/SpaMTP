.gene_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (!length(x)) return(character())
  x[is.na(x)] <- ""
  x <- sub("^(gene_symbol|symbol):", "symbol:", x)
  x <- sub("^(entrezgene|ncbigene|ncbi_gene):", "entrez:", x)
  x <- sub("^uniprotkb:", "uniprot:", x)
  x <- sub("^hgnc:hgnc:", "hgnc:", x)
  x[grepl("^[0-9]+$", x)] <- paste0("entrez:", x[grepl("^[0-9]+$", x)])
  ensembl <- grepl("^(ensembl:)?ensg[0-9]+(\\.[0-9]+)?$", x)
  x[ensembl] <- paste0("ensembl:", sub("\\.[0-9]+$", "", sub("^ensembl:", "", x[ensembl])))
  bare <- !grepl(":", x) & !grepl("^ramp_g_", x)
  x[bare] <- paste0("symbol:", x[bare])
  x
}

.gene_source_fingerprint <- function(source_df) {
  if (!is.data.frame(source_df) || !all(c("rampId", "sourceId") %in% names(source_df))) {
    stop("source_df must contain rampId and sourceId.", call. = FALSE)
  }
  keep <- !is.na(source_df$rampId) & grepl("^RAMP_G_", source_df$rampId)
  source <- data.frame(
    rampId = as.character(source_df$rampId[keep]),
    key = .gene_key(source_df$sourceId[keep]),
    commonName = if ("commonName" %in% names(source_df)) {
      tolower(trimws(as.character(source_df$commonName[keep])))
    } else rep(NA_character_, sum(keep))
  )
  source <- unique(source)
  source <- source[order(source$rampId, source$key, source$commonName,
                         na.last = TRUE, method = "radix"), , drop = FALSE]
  rownames(source) <- NULL
  # Hash identity-bearing contents, not source attributes or row ordering.
  # rlang is already an Imports dependency; no additional hash package is needed.
  rlang::hash(source)
}

.gene_validate_index <- function(index) {
  if (!inherits(index, "spamtp_gene_index")) {
    stop("index must be built by buildGeneMappingIndex().", call. = FALSE)
  }
  if (!identical(index$provenance$schema_version, 2L) ||
      !is.character(index$provenance$source_fingerprint) ||
      length(index$provenance$source_fingerprint) != 1L ||
      is.na(index$provenance$source_fingerprint) ||
      !nzchar(index$provenance$source_fingerprint) ||
      !is.character(index$conflict_gene_ids)) {
    stop("Gene index has an obsolete or incomplete schema. Rebuild it with buildGeneMappingIndex().",
         call. = FALSE)
  }
  invisible(TRUE)
}

.gene_check_source_index <- function(source_df, index) {
  .gene_validate_index(index)
  actual <- unique(as.character(source_df$rampId[
    !is.na(source_df$rampId) & grepl("^RAMP_G_", source_df$rampId)]))
  if (!setequal(actual, index$nodes$rampId) ||
      !identical(.gene_source_fingerprint(source_df), index$provenance$source_fingerprint)) {
    stop("gene_index does not match the RaMP source table identity fingerprint. ",
         "Rebuild it with buildGeneMappingIndex() using this source table.", call. = FALSE)
  }
  invisible(TRUE)
}

.gene_reference_lookup <- function(reference) {
  required <- c("hgnc_id", "symbol")
  if (!is.data.frame(reference) || !all(required %in% names(reference))) {
    stop("gene_reference must be an HGNC data frame with hgnc_id and symbol.", call. = FALSE)
  }
  reference <- as.data.frame(reference)
  if ("status" %in% names(reference)) {
    reference <- reference[!is.na(reference$status) & reference$status == "Approved", , drop = FALSE]
  }
  reference$hgnc_id <- toupper(trimws(reference$hgnc_id))
  if (!nrow(reference) || anyNA(reference$hgnc_id) || anyDuplicated(reference$hgnc_id) ||
      any(!grepl("^HGNC:[0-9]+$", reference$hgnc_id)) || anyNA(reference$symbol) ||
      any(!nzchar(trimws(reference$symbol))) || anyDuplicated(tolower(reference$symbol))) {
    stop("gene_reference must contain unique approved HGNC IDs and symbols.", call. = FALSE)
  }
  expand <- function(column, prefix, priority, method) {
    if (!column %in% names(reference)) return(NULL)
    values <- as.character(reference[[column]])
    values[is.na(values)] <- ""
    parts <- strsplit(values, "|", fixed = TRUE)
    ids <- rep(reference$hgnc_id, lengths(parts))
    values <- trimws(unlist(parts, use.names = FALSE))
    keep <- nzchar(values)
    if (!any(keep)) return(NULL)
    data.frame(key = .gene_key(paste0(prefix, values[keep])), gene_id = ids[keep],
               priority = rep(priority, sum(keep)), method = rep(method, sum(keep)))
  }
  lookup <- do.call(rbind, list(
    expand("hgnc_id", "", 0L, "hgnc_id"),
    expand("entrez_id", "entrez:", 0L, "entrez_id"),
    expand("ensembl_gene_id", "ensembl:", 0L, "ensembl_id"),
    expand("uniprot_ids", "uniprot:", 0L, "uniprot_id"),
    expand("symbol", "symbol:", 1L, "approved_symbol"),
    expand("prev_symbol", "symbol:", 2L, "previous_symbol"),
    expand("alias_symbol", "symbol:", 3L, "alias_symbol")
  ))
  # An approved symbol always wins over another gene's previous/alias symbol.
  lookup <- unique(lookup[order(lookup$key, lookup$priority, lookup$gene_id), ])
  best <- lookup$priority[match(lookup$key, lookup$key)]
  lookup <- lookup[lookup$priority == best, , drop = FALSE]
  list(reference = reference, lookup = lookup)
}

#' Build a gene identity index from RaMP and HGNC
#'
#' Connects RaMP gene records to HGNC using explicit stable identifiers, then
#' approved and previous symbols. Shared aliases are never used to merge RaMP
#' records. Conflicting stable identifiers quarantine the affected RaMP record.
#' Multiple non-conflicting RaMP records for one HGNC gene share one canonical
#' RaMP ID (lexicographically first); their pathway memberships can then be
#' united and counted once. No transitive merging across different HGNC genes
#' or orthology conversion is performed.
#' A fingerprint of normalized RaMP identifiers and common names guards index
#' reuse. Reordering rows, duplicating evidence or changing unrelated columns
#' does not invalidate it. Rebuild saved indices from older schema versions
#' with this function before using them in a new workflow.
#'
#' @param source_df RaMP identifier table from SpaMTPdb, with `sourceId`,
#'   `rampId` and optionally `commonName`.
#' @param gene_reference HGNC reference data frame. `NULL` retrieves the pinned
#'   reference through [SpaMTPdb::spaMTPdbGeneReference()]. Custom references
#'   require unique `hgnc_id` and `symbol`; optional fields are `status`,
#'   `entrez_id`, `ensembl_gene_id`, `uniprot_ids`, `prev_symbol`, `alias_symbol`.
#' @param gene_reference_version Version in SpaMTPdb's separate HGNC registry.
#' @param gene_reference_local_dir Local HGNC resource directory.
#' @param organism Must be `"Homo sapiens"`; this index is not a mouse-to-human
#'   orthology map. Use a separately curated species-specific database for other
#'   organisms.
#' @param offline Require local or cached HGNC data.
#' @return A `spamtp_gene_index` list containing reference lookup, gene groups,
#'   per-RaMP identity decisions, precomputed conflict membership, evidence and
#'   resource provenance including `source_fingerprint`.
#' @export
#' @examples
#' utils::str(formals(buildGeneMappingIndex))
buildGeneMappingIndex <- function(source_df, gene_reference = NULL,
                                  gene_reference_version = "latest",
                                  gene_reference_local_dir = NULL,
                                  organism = "Homo sapiens", offline = FALSE) {
  if (!identical(organism, "Homo sapiens")) {
    stop("HGNC gene mapping requires organism = 'Homo sapiens'; no orthology conversion is performed.",
         call. = FALSE)
  }
  if (!is.data.frame(source_df) || !all(c("rampId", "sourceId") %in% names(source_df))) {
    stop("source_df must contain rampId and sourceId.", call. = FALSE)
  }
  if (is.null(gene_reference)) {
    gene_reference <- SpaMTPdb::spaMTPdbGeneReference(
      version = gene_reference_version, local_dir = gene_reference_local_dir, offline = offline
    )
  }
  reference_provenance <- attr(gene_reference, "spamtp_gene_reference")
  if (!is.null(reference_provenance$organism) &&
      !identical(reference_provenance$organism, organism)) {
    stop("Gene reference organism does not match organism.", call. = FALSE)
  }
  prepared <- .gene_reference_lookup(gene_reference)
  reference <- prepared$reference
  lookup <- prepared$lookup
  source <- as.data.frame(source_df)
  source <- source[!is.na(source$rampId) & grepl("^RAMP_G_", source$rampId), , drop = FALSE]
  identifiers <- unique(data.frame(rampId = as.character(source$rampId),
                                   key = .gene_key(source$sourceId)))
  evidence <- merge(identifiers, lookup[lookup$priority < 3L, ], by = "key", sort = FALSE)
  if ("commonName" %in% names(source)) {
    labels <- unique(data.frame(rampId = as.character(source$rampId),
                                key = paste0("symbol:", tolower(trimws(source$commonName)))))
    extra <- merge(labels, lookup[lookup$priority %in% c(1L, 2L), ], by = "key", sort = FALSE)
    extra$priority <- extra$priority + 3L
    if (nrow(extra)) extra$method <- paste0("common_name_", extra$method)
    evidence <- rbind(evidence, extra)
  }
  evidence <- unique(evidence[order(evidence$rampId, evidence$priority, evidence$gene_id), ])
  top <- evidence[evidence$priority == evidence$priority[match(evidence$rampId, evidence$rampId)], ]
  choices <- lapply(split(seq_len(nrow(top)), top$rampId), function(rows) {
    identifiers <- lapply(split(top$gene_id[rows], top$key[rows]), unique)
    compatible <- Reduce(intersect, identifiers)
    # A shared protein accession can be disambiguated by a gene-specific
    # Entrez/Ensembl ID. Incompatible specific IDs still quarantine the record.
    if (length(compatible)) sort(compatible) else sort(unique(top$gene_id[rows]))
  })
  nodes <- data.frame(rampId = sort(unique(as.character(source$rampId))))
  candidates <- choices[nodes$rampId]
  nodes$candidate_gene_ids <- vapply(candidates, paste, character(1), collapse = ";")
  n <- lengths(candidates)
  nodes$gene_id <- ifelse(n == 1L, nodes$candidate_gene_ids, NA_character_)
  nodes$status <- ifelse(n == 0L, "unresolved_reference", ifelse(n == 1L, "resolved", "conflicting_record"))
  nodes$symbol <- reference$symbol[match(nodes$gene_id, reference$hgnc_id)]
  nodes$canonical_ramp_id <- nodes$rampId
  nodes$canonical_ramp_id[n > 1L] <- NA_character_
  groups <- split(nodes$rampId[n == 1L], nodes$gene_id[n == 1L])
  representatives <- vapply(groups, function(x) sort(x)[1L], character(1))
  nodes$canonical_ramp_id[n == 1L] <- unname(representatives[nodes$gene_id[n == 1L]])
  conflict_gene_ids <- sort(unique(as.character(unlist(strsplit(
    nodes$candidate_gene_ids[n > 1L], ";", fixed = TRUE), use.names = FALSE))))
  structure(list(nodes = nodes, lookup = lookup, reference = reference,
                 groups = groups, evidence = evidence, source_identifiers = identifiers,
                 conflict_gene_ids = conflict_gene_ids,
                 provenance = list(organism = organism, schema_version = 2L,
                                   source_fingerprint = .gene_source_fingerprint(source_df),
                                   source = attr(source_df, "spamtp_database"),
                                   reference = reference_provenance,
                                   custom_reference = is.null(reference_provenance))),
            class = "spamtp_gene_index")
}

.map_gene_index <- function(genes, index, ambiguous = c("exclude", "error"), verbose = TRUE) {
  ambiguous <- match.arg(ambiguous)
  .gene_validate_index(index)
  if ((!is.character(genes) && !is.numeric(genes)) || !is.null(dim(genes)) ||
      anyNA(genes) || any(!nzchar(trimws(as.character(genes))))) {
    stop("genes must be a vector without missing or blank identifiers.", call. = FALSE)
  }
  original_genes <- as.character(genes)
  original_keys <- .gene_key(original_genes)
  first <- !duplicated(original_keys)
  genes <- original_genes[first]
  keys <- original_keys[first]
  relevant <- index$lookup[index$lookup$key %in% keys, ]
  matches <- split(relevant, relevant$key)
  direct <- match(toupper(trimws(genes)), index$nodes$rampId)
  raw <- index$source_identifiers
  raw <- raw[raw$key %in% setdiff(keys, names(matches)), , drop = FALSE]
  raw_matches <- split(match(raw$rampId, index$nodes$rampId), raw$key)
  rows <- lapply(seq_along(genes), function(i) {
    gene_id <- symbol <- representative <- NA_character_
    ramp_ids <- candidates <- candidate_ramp_ids <- character()
    unresolved_competitors <- FALSE
    status <- "unrecognized_identifier"
    method <- NA_character_
    if (!is.na(direct[i])) {
      node <- index$nodes[direct[i], ]
      candidate_ramp_ids <- node$rampId
      method <- "ramp_id"
      if (node$status == "conflicting_record") {
        candidates <- strsplit(node$candidate_gene_ids, ";", fixed = TRUE)[[1L]]
        status <- "ambiguous"
      } else if (node$status == "unresolved_reference") {
        representative <- node$rampId
        ramp_ids <- node$rampId
        status <- "mapped_ramp_only"
      } else {
        candidates <- node$gene_id
      }
    } else if (!is.null(matches[[keys[i]]])) {
      hit <- matches[[keys[i]]]
      candidates <- sort(unique(hit$gene_id))
      method <- paste(unique(hit$method), collapse = ";")
    } else {
      # Preserve unique RaMP-only identifiers when HGNC has no entry, while
      # keeping their incomplete identity resolution visible in the audit.
      raw_nodes <- index$nodes[raw_matches[[keys[i]]], , drop = FALSE]
      resolved <- unique(raw_nodes$gene_id[!is.na(raw_nodes$gene_id)])
      candidate_ramp_ids <- sort(unique(raw_nodes$rampId))
      conflicts <- raw_nodes$candidate_gene_ids[raw_nodes$status == "conflicting_record"]
      candidates <- sort(unique(c(resolved, unlist(strsplit(conflicts, ";", fixed = TRUE),
                                                    use.names = FALSE))))
      if (nrow(raw_nodes)) method <- "ramp_crossreference"
      if (nrow(raw_nodes) == 1L && raw_nodes$status == "unresolved_reference") {
        representative <- raw_nodes$rampId
        ramp_ids <- representative
        status <- "mapped_ramp_only"
      } else if (nrow(raw_nodes) && any(raw_nodes$status != "resolved")) {
        # One resolved gene is not a unique match when competing records have
        # unknown or conflicting identities. Keep every raw candidate auditable.
        unresolved_competitors <- TRUE
        status <- "ambiguous"
      }
    }
    if (unresolved_competitors || length(candidates) > 1L) {
      status <- "ambiguous"
    } else if (length(candidates) == 1L) {
      gene_id <- candidates[1L]
      symbol <- index$reference$symbol[match(gene_id, index$reference$hgnc_id)]
      ramp_ids <- index$groups[[gene_id]]
      if (length(ramp_ids)) {
        representative <- sort(ramp_ids)[1L]
        status <- "mapped"
      } else {
        status <- "not_in_ramp"
        if (gene_id %in% index$conflict_gene_ids) status <- "conflicting_ramp_records"
      }
    }
    data.frame(input = genes[i], gene_id = gene_id, symbol = symbol,
               canonical_ramp_id = representative,
               ramp_ids = paste(sort(ramp_ids), collapse = ";"),
               status = status, method = method,
               candidate_gene_ids = paste(candidates, collapse = ";"),
               candidate_ramp_ids = paste(candidate_ramp_ids, collapse = ";"))
  })
  result <- if (length(rows)) do.call(rbind, rows) else data.frame(
    input = character(), gene_id = character(), symbol = character(),
    canonical_ramp_id = character(), ramp_ids = character(), status = character(),
    method = character(), candidate_gene_ids = character(), candidate_ramp_ids = character()
  )
  result <- result[match(original_keys, keys), , drop = FALSE]
  result$input <- original_genes
  rownames(result) <- NULL
  unresolved <- !result$status %in% c("mapped", "mapped_ramp_only")
  ambiguous_rows <- result$status %in% c("ambiguous", "conflicting_ramp_records")
  if (ambiguous == "error" && any(ambiguous_rows)) {
    stop("Ambiguous gene mapping: ", paste(unique(result$input[ambiguous_rows]), collapse = ", "), call. = FALSE)
  }
  if (verbose && any(unresolved)) warning(sum(unresolved), " gene input(s) unresolved; inspect mapping status and candidate_gene_ids.", call. = FALSE)
  attr(result, "gene_mapping") <- index$provenance
  result
}

#' Map gene identifiers with an auditable HGNC identity
#'
#' Accepts approved, previous and alias symbols (bare or `gene_symbol:`), HGNC,
#' Entrez, Ensembl gene IDs (optional version suffix), UniProt and RaMP IDs.
#' Ambiguous identifiers are excluded, never expanded to multiple genes.
#' Different RaMP records of the same resolved HGNC gene share one canonical
#' RaMP ID. Gene references are supplied by SpaMTPdb; this function does not
#' query live annotation services. Input order and duplicates are preserved.
#'
#' @param genes Vector of gene identifiers.
#' @param database Optional named SpaMTPdb resource bundle with `source_df`.
#' @param index Optional reusable index from [buildGeneMappingIndex()]. If
#'   provided, its recorded source and reference determine the mapping. If
#'   `database` is also supplied, its source-table fingerprint must match.
#' @param ambiguous Exclude ambiguous mappings with audit status, or error.
#' @param verbose Warn about unresolved input rows.
#' @param database_version RaMP version, independent of the HGNC version.
#' @param database_local_dir Local staged RaMP directory.
#' @inheritParams buildGeneMappingIndex
#' @return One row per input with `gene_id` (HGNC), approved `symbol`,
#'   `canonical_ramp_id`, all equivalent `ramp_ids` (semicolon separated),
#'   `status`, mapping `method` and `candidate_gene_ids`. The `gene_mapping`
#'   attribute records provenance. `mapped_ramp_only` explicitly denotes an
#'   unharmonized RaMP record without a resolved HGNC identity.
#'   `candidate_ramp_ids` records raw candidates encountered for direct RaMP
#'   IDs or source crossreferences, including unresolved/conflicting records.
#'   A shared crossreference with unresolved competitors remains `ambiguous`
#'   even if only one candidate has a resolved HGNC identity.
#' @export
#' @examples
#' utils::str(formals(mapGeneIdentifiers))
mapGeneIdentifiers <- function(genes, database = NULL, index = NULL,
                                gene_reference = NULL, gene_reference_version = "latest",
                                gene_reference_local_dir = NULL, organism = "Homo sapiens",
                                ambiguous = c("exclude", "error"), verbose = TRUE,
                                database_version = "latest", database_local_dir = NULL,
                                offline = FALSE) {
  if (is.null(index)) {
    source <- .spamtp_db_bundle("source_df", database = database,
                                version = database_version, local_dir = database_local_dir,
                                offline = offline)$source_df
    index <- buildGeneMappingIndex(source, gene_reference, gene_reference_version,
                                   gene_reference_local_dir, organism, offline)
  } else {
    .gene_validate_index(index)
    if (!identical(index$provenance$organism, organism)) {
      stop("Index organism does not match organism.", call. = FALSE)
    }
    if (!is.null(database)) {
      source <- .spamtp_db_bundle("source_df", database = database)$source_df
      .gene_check_source_index(source, index)
    }
  }
  .map_gene_index(genes, index, match.arg(ambiguous), verbose)
}

#' Store gene mapping in a native experiment's row metadata
#'
#' Adds aligned mapping columns to the selected experiment's `rowData()` and
#' provenance to its `metadata()`. Expression values, feature names, row order,
#' spatial coordinates and paired experiments are preserved. Multiple assay
#' rows that map to one gene are annotated, not silently summed; pathway
#' enrichment performs gene-level deduplication separately.
#'
#' @param object A SummarizedExperiment, SingleCellExperiment or
#'   SpatialExperiment, including native objects provided by SpaMTPData.
#' @param assay The primary experiment (`"main"`) or transcriptome altExp name.
#' @param id_column Gene identifier column in the selected `rowData()`.
#'   `NULL` uses its row names.
#' @param id_type Identifier type. `"auto"` accepts the same input formats as
#'   [mapGeneIdentifiers()]. Other values explicitly prefix unqualified IDs.
#' @param ... Passed to [mapGeneIdentifiers()], including `organism`, `index`,
#'   `database` and versioned reference settings.
#' @return The input object with `spamtp_gene_*` mapping columns and
#'   `metadata(experiment)$gene_mapping` provenance and input audit.
#' @export
#' @examples
#' utils::str(formals(annotateGeneIdentifiers))
annotateGeneIdentifiers <- function(object, assay = "main", id_column = NULL,
                                     id_type = c("auto", "symbol", "entrez", "ensembl", "uniprot", "hgnc"), ...) {
  id_type <- match.arg(id_type)
  args <- list(...)
  .gene_check_experiment_species(object, assay, args$organism %||% "Homo sapiens")
  experiment <- .experimentForAssay(object, assay)
  metadata <- SummarizedExperiment::rowData(experiment)
  if (is.null(id_column)) {
    ids <- rownames(experiment)
    if (is.null(ids)) stop("The selected experiment has no gene row names.", call. = FALSE)
  } else {
    if (length(id_column) != 1L || is.na(id_column) || !id_column %in% colnames(metadata)) {
      stop("id_column must name one column of the selected rowData().", call. = FALSE)
    }
    ids <- as.character(metadata[[id_column]])
  }
  if (id_type != "auto") {
    unqualified <- !is.na(ids) & !grepl(":", ids)
    ids[unqualified] <- paste0(id_type, ":", ids[unqualified])
  }
  mapping <- mapGeneIdentifiers(ids, ...)
  columns <- c(gene_id = "id", symbol = "symbol", canonical_ramp_id = "canonical_ramp_id",
                ramp_ids = "ramp_ids", status = "status", method = "method")
  for (column in names(columns)) metadata[[paste0("spamtp_gene_", columns[[column]])]] <- mapping[[column]]
  SummarizedExperiment::rowData(experiment) <- metadata
  S4Vectors::metadata(experiment)$gene_mapping <- list(
    provenance = attr(mapping, "gene_mapping"), id_column = id_column,
    id_type = id_type, inputs = mapping
  )
  .replaceExperiment(object, experiment, assay)
}

.gene_check_experiment_species <- function(object, assay, organism) {
  target <- .experimentForAssay(object, assay)
  metadata <- list(S4Vectors::metadata(object), S4Vectors::metadata(target))
  declared <- unique(unlist(lapply(metadata, function(x) c(x$SpaMTPData$organism, x$organism))))
  declared <- as.character(declared[!is.na(declared) & nzchar(declared)])
  if (length(declared) && any(declared != organism)) {
    stop("Experiment species (", paste(declared, collapse = ", "),
         ") does not match gene reference organism (", organism,
         "). No orthology conversion is performed.", call. = FALSE)
  }
  invisible(TRUE)
}

.gene_pathway_view <- function(resources, database = NULL, gene_mapping = c("auto", "hgnc", "ramp"),
                                gene_reference = NULL, gene_index = NULL,
                                gene_reference_version = "latest", gene_reference_local_dir = NULL,
                                organism = "Homo sapiens") {
  gene_mapping <- match.arg(gene_mapping)
  if (is.null(gene_reference) && !is.null(database$gene_reference)) gene_reference <- database$gene_reference
  official <- !is.null(attr(resources$source_df, "spamtp_database")) || is.null(database)
  enabled <- gene_mapping == "hgnc" || (gene_mapping == "auto" &&
    (official || !is.null(gene_reference) || !is.null(gene_index)))
  if (!enabled) return(list(resources = resources, index = NULL))
  if (is.null(gene_index)) gene_index <- buildGeneMappingIndex(
    resources$source_df, gene_reference, gene_reference_version, gene_reference_local_dir, organism
  )
  if (!inherits(gene_index, "spamtp_gene_index") || !identical(gene_index$provenance$organism, organism)) {
    stop("gene_index must be a matching organism index from buildGeneMappingIndex().", call. = FALSE)
  }
  .gene_check_source_index(resources$source_df, gene_index)
  canonicalize <- function(table) {
    table <- as.data.frame(table)
    i <- match(table$rampId, gene_index$nodes$rampId)
    hit <- !is.na(i)
    table$rampId <- as.character(table$rampId)
    table$rampId[hit] <- gene_index$nodes$canonical_ramp_id[i[hit]]
    table <- table[!is.na(table$rampId), , drop = FALSE]
    if ("commonName" %in% names(table)) {
      label <- gene_index$nodes$symbol[match(table$rampId, gene_index$nodes$rampId)]
      table$commonName <- as.character(table$commonName)
      table$commonName[!is.na(label)] <- label[!is.na(label)]
    }
    unique(table)
  }
  for (name in intersect(c("source_df", "analyte", "analytehaspathway"), names(resources))) {
    resources[[name]] <- canonicalize(resources[[name]])
  }
  list(resources = resources, index = gene_index)
}

.gene_prepare_de <- function(de, index, expand_records = FALSE) {
  audit <- .map_gene_index(as.character(de$gene), index)
  keep <- audit$status %in% c("mapped", "mapped_ramp_only")
  result <- de[keep, , drop = FALSE]
  audit_keep <- audit[keep, , drop = FALSE]
  result$rampId <- audit_keep$canonical_ramp_id
  result$commonName <- ifelse(is.na(audit_keep$symbol), result$gene, audit_keep$symbol)
  key <- paste(result$cluster, result$rampId, sep = "\r")
  duplicates <- split(seq_len(nrow(result)), key)
  for (rows in duplicates[lengths(duplicates) > 1L]) {
    values <- result[rows, intersect(c("avg_log2FC", "p_val_adj"), names(result)), drop = FALSE]
    if (nrow(unique(values)) > 1L) {
      stop("Multiple features map to gene ", result$commonName[rows[1L]],
           " with different differential statistics in cluster ", result$cluster[rows[1L]],
           ". Resolve or aggregate gene features before differential analysis.", call. = FALSE)
    }
  }
  first <- !duplicated(key)
  result <- result[first, , drop = FALSE]
  audit_keep <- audit_keep[first, , drop = FALSE]
  if (expand_records && nrow(result)) {
    ids <- strsplit(audit_keep$ramp_ids, ";", fixed = TRUE)
    result <- result[rep(seq_len(nrow(result)), lengths(ids)), , drop = FALSE]
    result$rampId <- unlist(ids, use.names = FALSE)
  }
  rownames(result) <- NULL
  attr(result, "gene_mapping") <- list(provenance = index$provenance, inputs = audit)
  result
}
