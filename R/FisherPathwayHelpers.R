.fisher_validate_input <- function(x, label) {
  allowed <- c("genes", "metabolites", "mzs")
  if (!is.list(x) || is.data.frame(x) || !length(x) ||
      is.null(names(x)) || anyNA(names(x)) || anyDuplicated(names(x)) ||
      !all(names(x) %in% allowed)) {
    stop(label, " must be a named list with unique keys: genes, metabolites, mzs.",
         call. = FALSE)
  }
  for (key in names(x)) {
    value <- x[[key]]
    if ((!is.null(value) && !is.character(value) && !is.numeric(value)) ||
        !is.null(dim(value)) || anyNA(value) ||
        any(!nzchar(trimws(as.character(value))))) {
      stop(label, "$", key, " must be a vector without missing or blank values.",
           call. = FALSE)
    }
    if (key == "mzs") {
      value <- suppressWarnings(as.numeric(sub(
        "^mz[-_]", "", trimws(as.character(value)), ignore.case = TRUE
      )))
      if (any(!is.finite(value) | value <= 0)) {
        stop(label, "$mzs must contain positive finite masses, optionally prefixed by mz- or mz_.",
             call. = FALSE)
      }
    } else {
      value <- trimws(as.character(value))
    }
    x[key] <- list(unique(value))
  }
  x
}

.fisher_modalities <- function(x) {
  unique(ifelse(names(x) == "genes", "G", "C"))
}

# Map source IDs first, then common names for otherwise unmatched queries.
# Keep all distinct RaMP mappings, but never count aliases as separate analytes.
.fisher_map_ids <- function(ids, type, source, known_ids) {
  query <- unique(tolower(ids))
  source <- source[grepl(paste0("^RAMP_", type, "_"), source$rampId), , drop = FALSE]
  key <- tolower(source$sourceId)
  selected <- key %in% query
  result <- data.frame(input = key[selected], rampId = source$rampId[selected])
  direct <- known_ids[grepl(paste0("^RAMP_", type, "_"), known_ids) &
                        tolower(known_ids) %in% query]
  result <- rbind(result, data.frame(input = tolower(direct), rampId = direct))
  if ("commonName" %in% names(source)) {
    key <- tolower(source$commonName)
    selected <- key %in% setdiff(query, result$input)
    result <- rbind(result, data.frame(input = key[selected], rampId = source$rampId[selected]))
  }
  unique(result)
}

.fisher_map_inputs <- function(foreground, universe, resources, args, verbose, gene_index = NULL) {
  source <- as.data.frame(resources$source_df)
  known_ids <- unique(c(source$rampId, resources$analyte$rampId,
                        resources$analytehaspathway$rampId))
  source <- source[!is.na(source$rampId) & nzchar(source$rampId), , drop = FALSE]
  known_ids <- known_ids[!is.na(known_ids) & nzchar(known_ids)]
  inputs <- list(Analyte = foreground, universe = universe)
  maps <- list()
  gene_audit <- NULL
  for (key in c("genes", "metabolites")) {
    ids <- unique(c(foreground[[key]], universe[[key]]))
    if (key == "genes" && !is.null(gene_index)) {
      gene_audit <- .map_gene_index(as.character(ids), gene_index, verbose = FALSE)
      keep <- gene_audit$status %in% c("mapped", "mapped_ramp_only")
      maps[[key]] <- data.frame(input = tolower(gene_audit$input[keep]),
                                rampId = gene_audit$canonical_ramp_id[keep])
    } else {
      maps[[key]] <- .fisher_map_ids(
        ids, if (key == "genes") "G" else "C", source, known_ids
      )
    }
  }
  masses <- unique(c(foreground$mzs, universe$mzs))
  mz_map <- data.frame(input = numeric(), rampId = character(), adduct = character())
  if (length(masses)) {
    verbose_message("Annotating foreground and universe m/z values together ...", verbose)
    if (is.null(args$db) && is.null(args$index)) args$db <- resources$chem_props
    annotations <- do.call(annotateTable, c(list(
      mz_df = data.frame(row_id = seq_along(masses), mz = masses),
      verbose = verbose
    ), args))
    # The current engine can also wrap legacy chemical databases without RaMP
    # IDs. Use source identifiers only for rows lacking a current RaMP ID.
    expanded <- .expand_annotation_column(annotations, "Ramp_IDs")
    expanded$Ramp_IDs <- toupper(expanded$Ramp_IDs)
    expanded <- expanded[expanded$Ramp_IDs %in% known_ids &
                           grepl("^RAMP_C_", expanded$Ramp_IDs), , drop = FALSE]
    mz_map <- data.frame(input = expanded$observed_mz, rampId = expanded$Ramp_IDs,
                         adduct = expanded$Adduct)
    legacy <- annotations[is.na(annotations$Ramp_IDs) |
                            !nzchar(trimws(annotations$Ramp_IDs)), , drop = FALSE]
    if (nrow(legacy)) {
      legacy <- .expand_annotation_column(legacy, "Isomers_IDs")
      mapping <- .fisher_map_ids(legacy$Isomers_IDs, "C", source, known_ids)
      legacy$input <- tolower(legacy$Isomers_IDs)
      legacy <- merge(legacy, mapping, by = "input")
      mz_map <- rbind(mz_map, data.frame(input = legacy$observed_mz,
                                        rampId = legacy$rampId, adduct = legacy$Adduct))
    }
  }
  mapped <- unmapped <- list()
  for (label in names(inputs)) {
    x <- inputs[[label]]
    rows <- data.frame(rampId = character(), sourceId = character(),
                       mz = numeric(), adduct = character())
    missing <- list()
    for (key in names(x)) {
      if (key == "mzs") {
        selected <- mz_map[mz_map$input %in% x[[key]], , drop = FALSE]
        rows <- rbind(rows, data.frame(rampId = selected$rampId,
                                      sourceId = selected$rampId,
                                      mz = selected$input, adduct = selected$adduct))
        missing[[key]] <- setdiff(x[[key]], selected$input)
      } else {
        selected <- maps[[key]][maps[[key]]$input %in% tolower(x[[key]]), , drop = FALSE]
        rows <- rbind(rows, data.frame(rampId = selected$rampId,
                                      sourceId = selected$input,
                                      mz = rep(NA_real_, nrow(selected)),
                                      adduct = rep(NA_character_, nrow(selected))))
        missing[[key]] <- x[[key]][!tolower(x[[key]]) %in% selected$input]
      }
      if (length(missing[[key]])) {
        warning(label, "$", key, ": ", length(missing[[key]]),
                " unmapped value(s) excluded; see attr(result, 'enrichment')$unmapped.",
                call. = FALSE)
      }
    }
    rows <- unique(rows)
    rows$commonName <- if ("commonName" %in% names(source)) {
      as.character(source$commonName[match(rows$rampId, source$rampId)])
    } else rep(NA_character_, nrow(rows))
    missing_name <- is.na(rows$commonName) | !nzchar(rows$commonName)
    rows$commonName[missing_name] <- rows$rampId[missing_name]
    mapped[[label]] <- rows
    unmapped[[label]] <- missing
  }
  list(mapped = mapped, unmapped = unmapped, gene_audit = gene_audit)
}

.fisher_pathway_details <- function(sets, foreground, mapped) {
  hits <- lapply(sets, intersect, y = foreground)
  collapse <- function(ids, type, column) {
    value <- mapped[mapped$rampId %in% ids &
                      grepl(paste0("^RAMP_", type, "_"), mapped$rampId), column]
    paste(unique(value[!is.na(value) & nzchar(value)]), collapse = ";")
  }
  result <- data.frame(row.names = seq_along(sets))
  for (type in c("C", "G")) {
    prefix <- if (type == "C") "metabolite" else "gene"
    result[[paste0(prefix, "_name_list")]] <- vapply(
      hits, collapse, character(1), type = type, column = "commonName"
    )
    result[[paste0(prefix, "_id_list")]] <- vapply(
      hits, collapse, character(1), type = type, column = "sourceId"
    )
  }
  result$adduct_info <- vapply(hits, function(ids) {
    x <- mapped[mapped$rampId %in% ids & !is.na(mapped$mz) &
                  !is.na(mapped$adduct), , drop = FALSE]
    if (!nrow(x)) return("")
    paste(unique(paste0(x$mz, "[", x$adduct, "]")), collapse = ";")
  }, character(1))
  result
}
