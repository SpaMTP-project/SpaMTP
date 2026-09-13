#' Fisher Exact Tests for Pathway Enrichment
#'
#' @param Analyte Foreground analytes as a named list with elements `genes`,
#'   `metabolites` and/or `mzs`. Identifiers are matched case-insensitively to
#'   source IDs (e.g. `gene_symbol:TP53`, `hmdb:HMDB0000122`), known RaMP IDs,
#'   then common names for otherwise unmatched queries. Duplicate mappings to
#'   the same RaMP ID count once. m/z values may be numeric or numeric strings,
#'   optionally prefixed with `mz-` or `mz_`.
#' @param max_path_size Maximum pathway size within the analysis universe
#'   (inclusive; default 500).
#' @param min_path_size Minimum pathway size within the analysis universe
#'   (inclusive; positive integer, default 5).
#' @param alternative Fisher test alternative: `"greater"`, `"less"` or
#'   `"two.sided"`.
#' @param pathway_all_info Include names, input IDs and m/z adduct information
#'   for foreground members of each pathway. Does not change the tests.
#' @param pval_cutoff Optional BH-adjusted p-value (FDR) cutoff in `[0, 1]`.
#'   Applied only after correcting over all eligible pathways.
#' @param verbose Whether to print progress messages. Mapping warnings are
#'   always emitted, including when `verbose = FALSE`.
#' @param database Optional named resource list from [loadSpaMTPDatabase()].
#'   Identifier analysis requires `source_df`, `analyte`, `analytehaspathway`
#'   and `pathway`; m/z analysis additionally uses `chem_props` unless `db`
#'   or `index` is supplied through `...`.
#' @param database_version SpaMTPdb/RaMP version used for pathway lookup.
#' @param database_source Database source; see [loadSpaMTPDatabase()].
#' @param database_local_dir Optional staged SpaMTPdb resource directory.
#' @param universe Optional measured background, in the same named-list format
#'   as `Analyte`. Supply all analytes eligible for foreground selection, not
#'   just significant hits. Foreground and universe must specify the same
#'   biological modalities (genes and/or compounds; `mzs` and `metabolites`
#'   both denote compounds). All mapped foreground IDs must belong to the
#'   mapped universe, otherwise an error is raised. Mapped background IDs
#'   without any pathway membership still contribute to the universe size.
#'   When `NULL`, use all pathway-linked IDs of the specified modalities in
#'   the database, before pathway-size filtering; mapped foreground IDs without
#'   pathway membership are then excluded with a warning.
#' @param gene_mapping `"auto"` uses HGNC identity harmonization for official
#'   SpaMTPdb resources and for custom bundles with `gene_reference` or
#'   `gene_index`. `"hgnc"` requires harmonization; `"ramp"` explicitly retains
#'   RaMP-only mapping for historical reproduction or a curated non-human
#'   database. Custom bundles without a reference use RaMP-only mapping in auto
#'   mode and never trigger a reference download.
#' @param gene_index Optional index built from the same source table by
#'   [buildGeneMappingIndex()]. Multiple RaMP records of a resolved HGNC gene
#'   have their pathway memberships united and count once. Conflicting RaMP
#'   records and ambiguous symbols are excluded and recorded in the audit.
#' @inheritParams buildGeneMappingIndex
#' @param ... Arguments passed to the indexed `annotateTable()` pipeline for
#'   m/z inputs, such as `ppm_error`, `adducts`, `db` or `index`.
#'
#' @details
#' With universe size U, foreground size K, pathway size M and overlap a, the
#' test table is `matrix(c(a, M - a, K - a, U - K - M + a), nrow = 2)`.
#' Memberships and identifiers are deduplicated by RaMP ID. Pathway sizes are
#' computed after intersecting with the universe; size filtering never changes
#' that universe. All pathways passing the size limits are tested, including
#' zero-overlap pathways. BH correction uses this complete family, before any
#' output cutoff. Missing display metadata does not remove a test.
#'
#' Unmapped input values are excluded with warnings and recorded in the
#' `enrichment` attribute. An empty mapped universe is an error. An empty
#' foreground gives p-values of 1 for all eligible pathways. If no pathways
#' pass the size limits, a zero-row data frame is returned.
#'
#' For a targeted panel, provide the measured panel as `universe`; the default
#' database background is not a substitute for the detection background.
#' m/z foreground and background values are annotated together with identical
#' settings. All mapped candidate RaMP IDs are retained; ambiguous mass matches
#' are not confirmed compound identifications. In a combined analysis, genes
#' and compounds are counted as individual RaMP IDs in a single pooled test.
#'
#' @return A data frame sorted by `p_val`, retaining the columns `pathway_name`,
#'   `pathway_id`, `type`, `pathwayCategory`, `p_val`, `fdr`, `ratio`,
#'   `analytes_in_pathways` (overlap) and `total_in_pathways` (size within the
#'   universe). Also includes `pathwayRampId`, `foreground_analytes_number`
#'   and `background_analytes_number`. `pathway_all_info = TRUE` adds member
#'   information. The `enrichment` attribute records mapped foreground and
#'   universe IDs, unmapped inputs, excluded foreground IDs, universe source,
#'   alternative, size limits and the number of tests before output filtering.
#'   Its `gene_mapping` entry contains per-input mapping status, gene groups,
#'   excluded conflicting RaMP records and versioned reference provenance.
#' @export
#' @import dplyr
#' @import stringr
#' @examples
#' utils::str(formals(fishersPathwayAnalysis))
#' # For a targeted gene panel:
#' # fishersPathwayAnalysis(
#' #   Analyte = list(genes = paste0("gene_symbol:", significant_genes)),
#' #   universe = list(genes = paste0("gene_symbol:", measured_panel_genes))
#' # )
fishersPathwayAnalysis <- function(Analyte,
                                  max_path_size = 500,
                                  min_path_size = 5,
                                  alternative = "greater",
                                  pathway_all_info = FALSE,
                                  pval_cutoff = NULL,
                                  verbose = TRUE,
                                  database = NULL,
                                  database_version = "latest",
                                  database_source = c("auto", "spamtpdb"),
                                  database_local_dir = NULL,
                                  universe = NULL,
                                  gene_mapping = c("auto", "hgnc", "ramp"),
                                  gene_reference = NULL,
                                  gene_index = NULL,
                                  gene_reference_version = "latest",
                                  gene_reference_local_dir = NULL,
                                  organism = "Homo sapiens",
                                  ...) {
  Analyte <- .fisher_validate_input(Analyte, "Analyte")
  if (!is.null(universe)) {
    universe <- .fisher_validate_input(universe, "universe")
    if (!setequal(.fisher_modalities(Analyte), .fisher_modalities(universe))) {
      stop("Analyte and universe must specify the same biological modalities (genes and/or compounds).",
           call. = FALSE)
    }
  }
  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  valid_size <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
  }
  if (!valid_size(min_path_size) || !valid_size(max_path_size) ||
      min_path_size > max_path_size) {
    stop("Pathway size limits must be positive finite integers with min_path_size <= max_path_size.",
         call. = FALSE)
  }
  if (!is.null(pval_cutoff) && (!is.numeric(pval_cutoff) ||
      length(pval_cutoff) != 1L || !is.finite(pval_cutoff) ||
      pval_cutoff < 0 || pval_cutoff > 1)) {
    stop("pval_cutoff must be NULL or one finite number in [0, 1].", call. = FALSE)
  }
  for (flag in c("verbose", "pathway_all_info")) {
    value <- get(flag)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(flag, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  args <- list(...)
  if (length(args) && (is.null(names(args)) || any(!nzchar(names(args))) ||
                       anyDuplicated(names(args)))) {
    stop("Annotation arguments in ... must have unique names.", call. = FALSE)
  }
  needed <- c("source_df", "analyte", "analytehaspathway", "pathway")
  if (length(c(Analyte$mzs, universe$mzs)) && is.null(args$db) && is.null(args$index)) {
    needed <- c(needed, "chem_props")
  }
  resources <- .spamtp_db_bundle(
    needed, database = database, version = database_version,
    source = match.arg(database_source), local_dir = database_local_dir
  )
  required <- list(source_df = c("sourceId", "rampId"), analyte = "rampId",
                   analytehaspathway = c("rampId", "pathwayRampId"),
                   pathway = "pathwayRampId")
  for (resource in names(required)) {
    if (!all(required[[resource]] %in% names(resources[[resource]]))) {
      stop(resource, " is missing required column(s): ",
           paste(setdiff(required[[resource]], names(resources[[resource]])), collapse = ", "),
           call. = FALSE)
    }
  }
  gene_mapping <- match.arg(gene_mapping)
  gene_view <- if ("genes" %in% names(Analyte)) .gene_pathway_view(
    resources, database, gene_mapping, gene_reference, gene_index,
    gene_reference_version, gene_reference_local_dir, organism
  ) else list(resources = resources, index = NULL)
  resources <- gene_view$resources
  mapping <- .fisher_map_inputs(Analyte, universe, resources, args, verbose, gene_view$index)
  foreground_ids <- unique(mapping$mapped$Analyte$rampId)
  links <- as.data.frame(resources$analytehaspathway)
  links <- unique(links[, c("rampId", "pathwayRampId"), drop = FALSE])
  pattern <- paste0("^RAMP_(", paste(.fisher_modalities(Analyte), collapse = "|"), ")_")
  links <- links[!is.na(links$rampId) & grepl(pattern, links$rampId) &
                   !is.na(links$pathwayRampId) & nzchar(links$pathwayRampId), , drop = FALSE]
  universe_ids <- if (is.null(universe)) unique(links$rampId) else {
    unique(mapping$mapped$universe$rampId)
  }
  if (!length(universe_ids)) stop("The mapped universe is empty.", call. = FALSE)
  outside <- setdiff(foreground_ids, universe_ids)
  if (length(outside)) {
    if (!is.null(universe)) {
      stop("Mapped foreground IDs are outside universe: ", paste(outside, collapse = ", "),
           ". Supply the complete measured background.", call. = FALSE)
    }
    warning(length(outside), " mapped foreground ID(s) without pathway membership excluded from the default universe.",
            call. = FALSE)
    foreground_ids <- intersect(foreground_ids, universe_ids)
  }
  links <- links[links$rampId %in% universe_ids, , drop = FALSE]
  sets <- split(as.character(links$rampId), as.character(links$pathwayRampId))
  sizes <- lengths(sets)
  sets <- sets[sizes >= min_path_size & sizes <= max_path_size]
  M <- lengths(sets)
  K <- length(foreground_ids)
  U <- length(universe_ids)
  a <- vapply(sets, function(ids) sum(ids %in% foreground_ids), integer(1))
  verbose_message(paste0("Testing ", length(sets), " pathways; foreground = ", K,
                          ", universe = ", U, "."), verbose)
  p <- vapply(seq_along(sets), function(i) {
    stats::fisher.test(matrix(c(a[i], M[i] - a[i], K - a[i],
                               U - K - M[i] + a[i]), nrow = 2),
                       alternative = alternative)$p.value
  }, numeric(1))
  metadata <- as.data.frame(resources$pathway)
  matched <- match(names(sets), metadata$pathwayRampId)
  field <- function(column, fallback = NA_character_) {
    value <- if (column %in% names(metadata)) as.character(metadata[[column]][matched]) else {
      rep(NA_character_, length(sets))
    }
    missing <- is.na(value) | !nzchar(value)
    value[missing] <- rep_len(fallback, length(sets))[missing]
    value
  }
  result <- data.frame(
    pathway_name = field("pathwayName", names(sets)),
    pathway_id = field("sourceId", names(sets)),
    type = field("type"), pathwayCategory = field("pathwayCategory"),
    p_val = p, fdr = stats::p.adjust(p, method = "BH"), ratio = unname(a / M),
    analytes_in_pathways = unname(a), total_in_pathways = unname(M),
    pathwayRampId = names(sets), foreground_analytes_number = rep(K, length(sets)),
    background_analytes_number = rep(U, length(sets)), stringsAsFactors = FALSE
  )
  if (pathway_all_info) {
    result <- cbind(result, .fisher_pathway_details(sets, foreground_ids, mapping$mapped$Analyte))
  }
  if (!is.null(pval_cutoff)) result <- result[result$fdr <= pval_cutoff, , drop = FALSE]
  result <- result[order(result$p_val, result$pathwayRampId), , drop = FALSE]
  rownames(result) <- NULL
  attr(result, "enrichment") <- list(
    universe_source = if (is.null(universe)) "database_pathway_members" else "measured",
    universe_ids = universe_ids, foreground_ids = foreground_ids,
    unmapped = mapping$unmapped, excluded_foreground_ids = outside,
    n_tested = length(sets), alternative = alternative,
    min_path_size = min_path_size, max_path_size = max_path_size,
    gene_mapping = if (is.null(gene_view$index)) list(mode = "ramp") else list(
      mode = "hgnc", provenance = gene_view$index$provenance, inputs = mapping$gene_audit,
      groups = gene_view$index$nodes[gene_view$index$nodes$canonical_ramp_id %in% universe_ids, ],
      conflicting_records = gene_view$index$nodes[gene_view$index$nodes$status == "conflicting_record", ]
    )
  )
  result
}

.expand_pathway_annotation_ids <- function(db_3) {
  if (!"Isomers_IDs" %in% names(db_3)) {
    stop("The stored legacy db_3 result must contain an Isomers_IDs column.")
  }
  db_3 <- tidyr::separate_rows(
    as.data.frame(db_3, stringsAsFactors = FALSE),
    dplyr::all_of("Isomers_IDs"),
    sep = "\\s*;\\s*"
  )
  db_3$Isomers_IDs <- trimws(as.character(db_3$Isomers_IDs))
  db_3[
    !is.na(db_3$Isomers_IDs) & nzchar(db_3$Isomers_IDs),
    ,
    drop = FALSE
  ]
}


#' Regional Pathway Enrichment
#'
#' This the function used to compute the gene/metabolites set enrichment for multi-omics spatial data
#'
#' @param SpaMTP A Bioconductor experiment contains spatial metabolomics(SM)/transcriptomics(ST) data or both, if contains SM data, it should be annotated via SpaMTP::annotateSM function.
#' @param ident A character scalar specifying the region column in cell metadata.
#' @param DE.list A list consisting of differential expression data.frames for each input modality. Within each data.frame column names MUST include 'cluster', 'gene', ('avg_log2FC' or 'logFC') and ('p_val_adj' or 'FDR').
#' @param analyte_types Vector of character strings defining which analyte types to use. Options can be c("genes"), c("metabolites") or both (default = c("genes", "metabolites")).
#' @param SM_assay Primary MSI experiment name (default = "main").
#' @param ST_assay Paired transcriptome altExp name (default = "transcriptome").
#' @param SM_slot The slot name containing the SM assay matrix data (default = "counts").
#' @param ST_slot The slot name containing the ST assay matrix data (default = "counts").
#' @param min_path_size The min number of metabolites in a specific pathway (default = 5).
#' @param max_path_size The max number of metabolites in a specific pathway (default = 500).
#' @param pval_cutoff_mets Adjusted p-value cutoff used when constructing
#'   metabolite ranks. Set to `1` to include annotated metabolites regardless
#'   of DE significance (default = 0.05).
#' @param pval_cutoff_genes A numerical value defining the adjusted p value cutoff for significant differentially expressed genes. If `NULL` cutoff = `0.05` (default = 0.05).
#' @param annotation_score_threshold Minimum indexed annotation score used for
#'   pathway mapping. It can be changed without re-running `annotateSM()` when
#'   the object was annotated with `min_score = 0` (default = 0.05).
#' @param annotation_source Metabolite annotation provenance. The default,
#'   `"current"`, requires scored RaMP IDs from the indexed annotation
#'   pipeline. `"auto"` permits a warned fallback to a legacy stored `db_3`,
#'   while `"legacy"` explicitly requests that compatibility path.
#' @param verbose Boolean indicating whether to show informative messages. If FALSE these messages will be suppressed (default = TRUE).
#' @param database Optional named list of database resources, normally created
#'   by [loadSpaMTPDatabase()].
#' @param database_version SpaMTPdb/RaMP version used for pathway lookup.
#' @param database_source Database source; see [loadSpaMTPDatabase()].
#' @param database_local_dir Optional staged SpaMTPdb resource directory.
#'
#' @return A SpaMTP object with set enrichment on given analyte types.
#' @export
#'
#' @importFrom rlang %||%
#'
#' @examples
#' utils::str(formals(findRegionalPathways))
#' # SpaMTP = findRegionalPathways(SpaMTP, polarity = "positive")
#' @inheritParams fishersPathwayAnalysis
#' @details Gene mapping uses the same identity index as
#'   [fishersPathwayAnalysis()]. Multiple RaMP records of a gene have their
#'   memberships united. If multiple input features map to that gene with
#'   different differential statistics in one cluster, resolve the duplicate
#'   features before differential analysis; the function does not select an
#'   arbitrary statistic.
findRegionalPathways = function(SpaMTP,
                                ident,
                                DE.list,
                                analyte_types = c("genes", "metabolites"),
                                SM_assay = "main",
                                ST_assay = "transcriptome",
                                SM_slot = "counts",
                                ST_slot = "counts",
                                min_path_size = 5,
                                max_path_size = 500,
                                pval_cutoff_mets = 0.05,
                                pval_cutoff_genes = 0.05,
                                annotation_score_threshold = 0.05,
                                annotation_source = c("current", "auto", "legacy"),
                                verbose = TRUE,
                                database = NULL,
                                database_version = "latest",
                                database_source = c("auto", "spamtpdb"),
                                database_local_dir = NULL,
                                gene_mapping = c("auto", "hgnc", "ramp"),
                                gene_reference = NULL, gene_index = NULL,
                                gene_reference_version = "latest",
                                gene_reference_local_dir = NULL,
                                organism = "Homo sapiens") {
  annotation_source <- match.arg(annotation_source)
  database_resources <- .spamtp_db_bundle(
    c("chem_props", "source_df", "analytehaspathway", "pathway"),
    database = database,
    version = database_version,
    source = match.arg(database_source),
    local_dir = database_local_dir
  )
  gene_view <- if ("genes" %in% analyte_types) .gene_pathway_view(
    database_resources, database, gene_mapping, gene_reference, gene_index,
    gene_reference_version, gene_reference_local_dir, organism
  ) else list(resources = database_resources, index = NULL)
  database_resources <- gene_view$resources
  if (!is.null(gene_view$index)) .gene_check_experiment_species(SpaMTP, ST_assay, organism)
  chem_props <- database_resources$chem_props
  source_df <- database_resources$source_df
  analytehaspathway <- database_resources$analytehaspathway
  pathway <- database_resources$pathway

  ## Checks for ident in SpaMTP Object
  cellMetadata <- .cellMetadata(SpaMTP)
  if (!(ident %in% colnames(cellMetadata))) {
    stop(
      "Ident: ",
      ident,
      " not found in the SpaMTP object's cell metadata. Make sure the column exists and is a factor."
    )
  }
  cluster_vector = as.factor(cellMetadata[[ident]])
  assignment = cluster_vector
  cluster = levels(cluster_vector)
  ## Checks for data in SM and/or ST assay
  if ("genes" %in% analyte_types) {
    st_obj <- .assayData(SpaMTP, ST_assay, ST_slot)
    if (is.null(st_obj)) {
      stop(
        paste0(
          "No data exists in object[[",
          ST_assay,
          "]][",
          ST_slot,
          "] .. If you are using transcriptomic data with 'genes' in 'analyte_types', please ensure this dataslot exists within your SpaMTP object, else remove 'genes' from analyte_tpes"
        )
      )
    } else{
      gene_matrix = Matrix::t(st_obj)
      if (length(cluster_vector) != nrow(gene_matrix)) {
        stop(
          "Please make sure the input ident is a vector the same length as the number of spots/cells in the gene assay!"
        )
      }
    }
  }
  if ("metabolites" %in% analyte_types) {
    sm_obj <- .assayData(SpaMTP, SM_assay, SM_slot)
    if (is.null(sm_obj)) {
      stop(
        paste0(
          "No data exists in object[[",
          SM_assay,
          "]][",
          SM_slot,
          "] .. If you are using metabolic data with 'metabolites' in 'analyte_types', please ensure this dataslot exists within your SpaMTP object, else remove 'metabolites' from analyte_tpes"
        )
      )
    } else{
      mass_matrix = Matrix::t(sm_obj)
      if (length(cluster_vector) != nrow(mass_matrix)) {
        stop(
          "Please make sure the input ident is a vector the same length as the number of spots/cells in the metabolite assay!"
        )
      }
    }
  }
  annotation_metadata <- NULL
  db_3 <- NULL
  if ("metabolites" %in% analyte_types) {
    # (2) Resolve annotations. Scored Ramp_IDs from the indexed pipeline are
    # the primary key; source-ID joins are only available in legacy mode.
    verbose_message(
      message_text = "Resolving current RaMP metabolite annotations ... ",
      verbose = verbose
    )
    db_3 <- .resolve_pathway_metabolite_annotations(
      SpaMTP,
      annotation_source = annotation_source,
      score_threshold = annotation_score_threshold,
      chemical_properties = chem_props
    )
    annotation_metadata <- attr(db_3, "annotation_metadata")
    annotation_metadata$pathway_pval_cutoff_mets <- pval_cutoff_mets
    if (isTRUE(verbose) && length(annotation_metadata$engine)) {
      ramp_label <- annotation_metadata$ramp_version
      if (is.null(ramp_label) || is.na(ramp_label)) ramp_label <- "unknown"
      message(
        "Annotation engine: ", annotation_metadata$engine,
        "; RaMP: ", ramp_label
      )
    }
  }

  ### Adding DE Results
  if (length(DE.list) != length(analyte_types)) {
    stop(
      "Number of DE data.frames provided does not match the number of analyte types specified. Please make sure a DE dataframe is provided for each analyte type"
    )
  }
  verbose_message(message_text = "Constructing DE dataframes.... ", verbose = verbose)
  for (i in 1:length(analyte_types)) {
    verbose_message(
      message_text = paste0(
        "Assuming DE.list[",
        i,
        "] contains ",
        analyte_types[i] ,
        " results .... "
      ),
      verbose = verbose
    )
    DE <- DE.list[[i]]
    if (any(c("avg_log2FC", "logFC") %in% colnames(DE)) &&
        any(c("p_val_adj", "FDR") %in% colnames(DE)) &&
        "cluster" %in% colnames(DE) &&
        "gene" %in% colnames(DE)) {
      if ("logFC" %in% colnames(DE)) {
        colnames(DE)[colnames(DE) == "logFC"] <- "avg_log2FC"
      }
      # Rename FDR to p_val_adj if FDR exists
      if ("FDR" %in% colnames(DE)) {
        colnames(DE)[colnames(DE) == "FDR"] <- "p_val_adj"
      }
      if (analyte_types[i] == "metabolites") {
        DE = DE %>% rename(mz_name = gene)
        db_3 = merge(db_3 , DE, by = "mz_name")
        DE.list[[analyte_types[i]]] <- db_3
      } else {
        if (!is.null(gene_view$index)) {
          source_gene <- .gene_prepare_de(DE, gene_view$index)
        } else {
          DE = DE %>% mutate(commonName = toupper(gene))
          source_gene = merge(DE, source_df[which(grepl(source_df$rampId, pattern = "RAMP_G")), ], by = "commonName")
        }
        DE.list[[analyte_types[i]]] <- source_gene
      }
    } else {
      stop(
        "DE dataframe [",
        i,
        "] provided does not have the correct column names ... column names MUST include 'cluster', 'gene', ('avg_log2FC' or 'logFC') and ('p_val_adj' or 'FDR'). Please adjust column names in all DE data.frames to match ..."
      )
    }
  }
  # Get pathway db
  verbose_message(message_text = "Constructing pathway database ..." , verbose = verbose)
  chempathway = merge(analytehaspathway, pathway, by = "pathwayRampId")

  pathway_db = lapply(split(chempathway$rampId, chempathway$pathwayName), unique)
  pathway_db = pathway_db[which(!duplicated(tolower(names(pathway_db))))]
  pathway_db = pathway_db[lapply(pathway_db, length) >= min_path_size  &
                            lapply(pathway_db, length) <= max_path_size]

  gc()
  gsea_all_cluster = data.frame()
  all_ranks = list()
  pb3 = txtProgressBar(
    min = 0,
    max = length(cluster),
    initial = 0,
    style = 3
  )
  for (i in cluster) {
    i <- as.character(i)
    ranks <- c()
    if ("metabolites" %in% analyte_types) {
      ## metabolites
      DE_met <- DE.list[["metabolites"]]
      sub_db3 = DE_met[which(as.character(DE_met$cluster) == i), ] %>% dplyr::filter(p_val_adj <= pval_cutoff_mets %||% 0.05) %>% dplyr::filter(!duplicated(ramp_id))
      met_ranks = scale(sub_db3$avg_log2FC, center = 0)
      names(met_ranks) = sub_db3$ramp_id
      ranks <- c(ranks, met_ranks)
    }
    if ("genes" %in% names(DE.list)) {
      ## genes
      DE_rna <- DE.list[["genes"]]
      sub_de_gene = DE_rna[which(as.character(DE_rna$cluster) == i), ] %>% dplyr::filter(p_val_adj <= pval_cutoff_genes %||% 0.05) %>% dplyr::filter(!duplicated(rampId))
      ranks_gene_vector = scale(sub_de_gene$avg_log2FC, center = 0)
      names(ranks_gene_vector) = sub_de_gene$rampId
      # Genes and metabolites
      ranks <- c(ranks, ranks_gene_vector)
    }

    ranks = ranks[which(!duplicated(names(ranks)))]
    all_ranks[[i]] = ranks[is.finite(ranks)]

    gsea_result <- c()
    if (length(all_ranks[[i]]) > 0) {
      suppressWarnings({
        gsea_result = fgsea::fgsea(
          pathways =  pathway_db,
          stats = all_ranks[[i]],
          minSize = min_path_size,
          maxSize = max_path_size
        )  %>%  dplyr::mutate(Cluster_id = i)
      })

    } else {
      gsea_result <- data.table::data.table(
        pathway = character(0),
        pval = numeric(0),
        padj = numeric(0),
        log2err = numeric(0),
        ES = numeric(0),
        NES = numeric(0),
        size = integer(0),
        leadingEdge = list(),
        Cluster_id = i
      )
    }
    gsea_result = na.omit(gsea_result) %>% filter(!duplicated(pathway))
    short_source = source_df[which((source_df$rampId %in% names(all_ranks[[i]])) &
                                     !duplicated(source_df$rampId)), ]

    if (!nrow(gsea_result)) next
    addtional_entry = do.call(rbind, lapply(seq_len(nrow(gsea_result)), function(x) {
      temp = unique(unlist(gsea_result$leadingEdge[x]))
      if ("metabolites" %in% analyte_types) {
        temp_ref =   sub_db3[which(sub_db3$ramp_id %in% temp), ] %>% dplyr::mutate(adduct_info = paste0(observed_mz, "[", Adduct, "]")) %>% dplyr::filter(!duplicated(adduct_info))
      }
      if ("genes" %in% analyte_types) {
        temp_rna = short_source[which((short_source$rampId %in% temp) &
                                        (grepl(short_source$rampId, pattern = "RAMP_G"))), ]
      }
      return(
        data.frame(
          adduct_info = if("metabolites" %in% analyte_types){paste0(temp_ref$adduct_info, collapse = ";")}else{""},
          leadingEdge_metabolites = if("metabolites" %in% analyte_types){paste0(sub(";.*", "", temp_ref$IsomerNames), collapse = ";")}else{""},
          leadingEdge_metabolites_id = if("metabolites" %in% analyte_types){paste0(temp_ref$chem_source_id, collapse = ";")}else{""},
          leadingEdge_genes = if("genes" %in% analyte_types){paste0(temp_rna$commonName, collapse = ";")}else{""},
          met_regulation = if("metabolites" %in% analyte_types){paste0(ifelse(ranks[which((names(ranks) %in% temp) &
                                                       (grepl(names(ranks), pattern = "RAMP_C")))] >= 0, "\u2191", "\u2193"), collapse = ";")}else{""},
          rna_regulation = if("genes" %in% names(DE.list)){paste0(ifelse(ranks[which((names(ranks) %in% temp) &
                                                       (grepl(names(ranks), pattern = "RAMP_G")))] >= 0, "\u2191", "\u2193"), collapse = ";")}else{""}
        )
      )
    }))
    gsea_result = cbind(gsea_result , addtional_entry)
    gsea_all_cluster = rbind(gsea_all_cluster, gsea_result)
    setTxtProgressBar(pb3, as.numeric(which(cluster == i)))
  }
  close(pb3)

  if (!nrow(gsea_all_cluster)) {
    result <- data.frame(pathwayName = character(), pval = numeric(), padj = numeric(),
                          NES = numeric(), Cluster_id = character())
    attr(result, "annotation_metadata") <- annotation_metadata
    attr(result, "gene_mapping") <- if (!is.null(gene_view$index)) attr(DE.list[["genes"]], "gene_mapping") else list(mode = "ramp")
    return(result)
  }
  gsea_all_cluster <- na.omit(gsea_all_cluster)%>%
    dplyr::mutate(group_importance = sum(abs(NES)))
  colnames(gsea_all_cluster)[1] = "pathwayName"
  gsea_all_cluster = merge(gsea_all_cluster, pathway, by = "pathwayName")
  attr(gsea_all_cluster, "annotation_metadata") <- annotation_metadata
  attr(gsea_all_cluster, "gene_mapping") <- if (!is.null(gene_view$index)) {
    attr(DE.list[["genes"]], "gene_mapping")
  } else list(mode = "ramp")
  return(gsea_all_cluster)
}




#' Runs multilevel Monte-Carlo variant for performing gene sets co-regulation analysis using the RAMP_DB metabolite/gene database.
#'
#' This function is adapted from the [fgsea::geseca](https://github.com/alserglab/fgsea/blob/master/R/geseca-multilevel.R) package to identify significantly expressed RAMP_DB pathways based on an expression/feature embedding matrix.
#'
#' @param E expression matrix, rows corresponds to RAMP_IDs, columns corresponds to cell barcodes.
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded (default = 1).
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded (default = `nrow(E) - 1`).
#' @param center a logical value indicating whether the gene expression should be centered to have zero mean before the analysis takes place (default = TRUE).
#' @param scale a logical value indicating whether the gene expression should be scaled to have unit variance before the analysis takes place (default = FALSE).
#' @param sampleSize sample size for conditional sampling (default = 101).
#' @param eps This parameter sets the boundary for calculating P-values (default = 1e-50).
#' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 0).
#' @param BPPARAM Parallelization parameter used in bplapply (default = NULL).
#' @param nPermSimple Number of permutations in the simple geseca implementation for preliminary estimation of P-values (default = 1000).
#' @param database Optional named list of database resources, normally created
#'   by [loadSpaMTPDatabase()].
#' @param database_version SpaMTPdb/RaMP version used for pathway lookup.
#' @param database_source Database source; see [loadSpaMTPDatabase()].
#' @param database_local_dir Optional staged SpaMTPdb resource directory.
#'
#' @return A table with GESECA results. Each row corresponds to a tested RAMP_DB pathway.
#' @export
#'
#' @examples
#' utils::str(formals(runRAMPGeseca))
#' # sig_pathways <- runRAMPGeseca(E, minSize=15, maxSize=500)
runRAMPGeseca <- function(E,
                          minSize     = 1,
                          maxSize     = nrow(E) - 1,
                          center      = TRUE,
                          scale       = FALSE,
                          sampleSize  = 101,
                          eps         = 1e-50,
                          nproc       = 0,
                          BPPARAM     = NULL,
                          nPermSimple = 1000,
                          database = NULL,
                          database_version = "latest",
                          database_source = c("auto", "spamtpdb"),
                          database_local_dir = NULL){

  database_resources <- .spamtp_db_bundle(
    c("analytehaspathway", "pathway"),
    database = database,
    version = database_version,
    source = match.arg(database_source),
    local_dir = database_local_dir
  )
  analytehaspathway <- database_resources$analytehaspathway
  pathway <- database_resources$pathway

  chempathway = merge(analytehaspathway, pathway, by = "pathwayRampId")

  pathway_db = split(chempathway$rampId, chempathway$pathwayName)
  pathway_db = pathway_db[which(!duplicated(tolower(names(pathway_db))))]
  pathway_db = pathway_db[lapply(pathway_db, length) >= minSize  &
                            lapply(pathway_db, length) <= maxSize]

  gesecaRes <- fgsea::geseca(pathway_db, E, minSize = minSize, maxSize = maxSize, center = center, scale = scale,sampleSize = sampleSize, eps = eps, nproc = nproc, BPPARAM = BPPARAM, nPermSimple = nPermSimple)

  return(gesecaRes)

}


#' Create a Pathway Assay from Gene or Metabolite Data
#'
#' This function creates a new assay within the provided Bioconductor experiment which contains features (either genes or metabolites) labeled by their respective RAMP ID. This assay can be used for running feature set co-regulation analysis (based on GSCA; \doi{10.1093/bioinformatics/btp502}).
#'
#' @param SpaMTP A Bioconductor experiment containing either spatial metabolic or transcriptomic data
#' @param analyte_type Character string specifying the type of analytes to process.Must be either "genes" or "metabolites" (default = "metabolites").
#' @param assay Character string specifying the name of the assay to use as source data (default = "Spatial").
#' @param slot Character string specifying which slot in the assay to use as source data (default = "counts").
#' @param new_assay Character string specifying the name of the new assay to create (default = "pathway").
#' @param annotation_score_threshold Minimum indexed annotation score used to
#'   map m/z features to RaMP compounds (default = 0.05).
#' @param annotation_source Metabolite annotation provenance. `"current"`
#'   requires the indexed, scored RaMP output; `"auto"` and `"legacy"` enable
#'   compatibility with older serialized SpaMTP objects.
#' @param verbose Boolean logical value indicating whether to print verbose messages during execution. (default = TRUE).
#' @param database Optional named list of database resources, normally created
#'   by [loadSpaMTPDatabase()].
#' @param database_version SpaMTPdb/RaMP version used for pathway lookup.
#' @param database_source Database source; see [loadSpaMTPDatabase()].
#' @param database_local_dir Optional staged SpaMTPdb resource directory.
#'
#' @return A SpaMTP object with a new assay added, containing respective gene/metabolite data formatted based on RAMP_db IDs.
#' @export
#'
#' @importFrom dplyr mutate group_by summarise ungroup
#' @importFrom tidyr separate_rows
#' @importFrom data.table as.data.table
#'
#' @examples
#' utils::str(formals(createPathwayAssay))
#' ## Create a pathway assay from metabolite data
#' #spamtp_obj <- createPathwayAssay(spamtp_obj, analyte_type = "metabolites", assay = "SPM", new_assay = "pathway")
#'
#' ## Create a pathway assay from gene data with verbose output
#' #spamtp_obj <- createPathwayAssay(spamtp_obj, analyte_type = "genes", assay = "SPT", new_assay = "gene_pathway", verbose = TRUE)
createPathwayAssay <- function(SpaMTP, analyte_type = "metabolites", assay = "Spatial", slot = "counts", new_assay = "pathway", annotation_score_threshold = 0.05, annotation_source = c("current", "auto", "legacy"), verbose = TRUE, database = NULL, database_version = "latest", database_source = c("auto", "spamtpdb"), database_local_dir = NULL){
  .requireExperiment(SpaMTP, "SingleCellExperiment")


  annotation_source <- match.arg(annotation_source)
  database_resources <- .spamtp_db_bundle(
    c("chem_props", "source_df"),
    database = database,
    version = database_version,
    source = match.arg(database_source),
    local_dir = database_local_dir
  )
  chem_props <- database_resources$chem_props
  source_df <- database_resources$source_df

  if(!analyte_type %in% c("genes", "metabolites")){
    stop("Incorrect `analyte_type` provided! must be either 'genes' or 'metabolites'. Please provided the correct analyte matching the selected assay data.")
  }

  if (analyte_type == "genes") {
    assayMatrix <- tryCatch(.assayData(SpaMTP, assay, slot), error = function(e) NULL)
    if (is.null(assayMatrix)) {
      stop(
        paste0(
          "No data exists in object[[",
          assay,
          "]][",
          slot,
          "] .. If you are using transcriptomic data with 'genes' in 'analyte_types', please ensure this dataslot exists within your SpaMTP object, else remove 'genes' from analyte_tpes"
        )
      )
    } else{
      matrix <- as.data.frame(assayMatrix)
      matrix$commonName <- toupper(rownames(matrix))
      matrix = merge(matrix,
                     unique(source_df[which(grepl(source_df$rampId, pattern = "RAMP_G")), ][c("rampId" ,"commonName")]),
                     by = "commonName")
      dupe_list <- split(which(matrix$rampId %in% matrix$rampId[duplicated(matrix$rampId)]),
                         matrix$rampId[matrix$rampId %in% matrix$rampId[duplicated(matrix$rampId)]])

      meta.data <- matrix[c("rampId" ,"commonName")] %>%
        group_by(rampId) %>%
        summarise(commonName = paste(commonName, collapse = "; ")) %>%
        ungroup()


      matrix$commonName <- NULL

    }
  }
  if (analyte_type == "metabolites") {
    assayMatrix <- tryCatch(.assayData(SpaMTP, assay, slot), error = function(e) NULL)
    if (is.null(assayMatrix)) {
      stop(
        paste0(
          "No data exists in object[[",
          assay,
          "]][",
          slot,
          "] .. If you are using metabolic data with 'metabolites' in 'analyte_types', please ensure this dataslot exists within your SpaMTP object, else remove 'metabolites' from analyte_tpes"
        )
      )
    } else{
      matrix <- as.data.frame(assayMatrix)

      # (2) Annotation
      verbose_message(
        message_text = "Resolving current RaMP metabolite annotations ... ",
        verbose = verbose
      )
      db_3 <- .resolve_pathway_metabolite_annotations(
        SpaMTP,
        annotation_source = annotation_source,
        score_threshold = annotation_score_threshold,
        chemical_properties = chem_props
      )

      ### Adding DE Results
      db_3 <- db_3[c("mz_name",  "ramp_id")]
      db_3 <- db_3 %>% distinct()
      matrix$mz_name <- rownames(assayMatrix)
      matrix = merge(db_3 , matrix, by = "mz_name")

      meta.data <- matrix[c("ramp_id" ,"mz_name")] %>%
        group_by(ramp_id) %>%
        summarise(mz_name = paste(mz_name, collapse = "; ")) %>%
        ungroup()

      meta.data$rampId <- meta.data$ramp_id
      meta.data <- meta.data[c("rampId","mz_name")]

      rm(db_3)

      dupe_list <- split(which(matrix$ramp_id %in% matrix$ramp_id[duplicated(matrix$ramp_id)]),
                         matrix$ramp_id[matrix$ramp_id %in% matrix$ramp_id[duplicated(matrix$ramp_id)]])

      matrix$mz_name <- NULL
      matrix$rampId <- matrix$ramp_id
      matrix$ramp_id <- NULL

    }
  }


  # Aggregate by explicit RaMP IDs, never by data.table row names. Select
  # expression columns in their original order after database joins.
  identifiers <- unique(matrix$rampId)
  if (!length(identifiers)) {
    stop("No features map to RaMP IDs in the selected database.", call. = FALSE)
  }
  values <- as.matrix(matrix[, colnames(assayMatrix), drop = FALSE])
  storage.mode(values) <- "numeric"
  if (any(!is.finite(values))) {
    stop("RaMP aggregation requires finite expression values.", call. = FALSE)
  }
  group <- match(matrix$rampId, identifiers)
  weights <- Matrix::sparseMatrix(
    i = group, j = seq_along(group),
    x = 1 / tabulate(group, length(identifiers))[group],
    dims = c(length(identifiers), nrow(values)))
  pathwayMatrix <- weights %*% values
  dimnames(pathwayMatrix) <- list(identifiers, colnames(assayMatrix))
  pathwayMetadata <- data.frame(
    rampId = rownames(pathwayMatrix),
    row.names = rownames(pathwayMatrix)
  )
  pathwayMetadata <- merge(pathwayMetadata, meta.data, by = "rampId", all = TRUE)
  rownames(pathwayMetadata) <- pathwayMetadata$rampId
  pathwayMetadata <- pathwayMetadata[rownames(pathwayMatrix), , drop = FALSE]

  pathwayExperiment <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = pathwayMatrix),
    rowData = S4Vectors::DataFrame(pathwayMetadata)
  )
  SingleCellExperiment::altExp(SpaMTP, new_assay) <- pathwayExperiment


  return(SpaMTP)

}



############################# Annotation Estimation using Pathway Results ########################################



#' Store pathway scores in an alternative experiment
#'
#' This function computes pathway-level scores from analyte-level expression data
#' and stores them in altExp for SingleCellExperiment-derived objects.
#' Each score is the sum of matched expression divided by the square root
#' of the matched feature count, then centred and scaled across pixels.
#' Constant scores are zero; unmatched pathways are removed or stored as NA.
#'
#' @param object A SingleCellExperiment (including SpatialExperiment) with expression
#'   features indexed by RaMP IDs.
#' @param assay Primary or alternative experiment to score.
#' @param slot Expression assay. If omitted for a Bioconductor container,
#'   prefers logcounts, normcounts, then counts in the selected experiment.
#' @param new.assay Character. Name of the new assay where pathway scores will be stored (defaults = "pathway").
#' @param remove.nans Remove pathways without matched analytes (default TRUE).
#' @param database Optional named list of database resources, normally created
#'   by [loadSpaMTPDatabase()].
#' @param database_version SpaMTPdb/RaMP version used for pathway lookup.
#' @param database_source Database source; see [loadSpaMTPDatabase()].
#' @param database_local_dir Optional staged SpaMTPdb resource directory.
#'
#' @return The input with pathwayScores in a new alternative experiment.
#'   Pathway identifiers are preserved.
#'
#' @export
#'
#' @examples
#' utils::str(formals(createPathwayObject))
createPathwayObject <- function(object,
                                assay = NULL,
                                slot = "logcounts",
                                new.assay = "pathway",
                                remove.nans = TRUE,
                                database = NULL,
                                database_version = "latest",
                                database_source = c("auto", "spamtpdb"),
                                database_local_dir = NULL
) {
  .requireExperiment(object, "SingleCellExperiment")


  database_resources <- .spamtp_db_bundle(
    c("analytehaspathway", "pathway"),
    database = database,
    version = database_version,
    source = match.arg(database_source),
    local_dir = database_local_dir
  )
  analytehaspathway <- database_resources$analytehaspathway
  pathway <- database_resources$pathway

  chempathway = merge(analytehaspathway, pathway, by = "pathwayRampId")
  pathway_db <- split(chempathway$rampId, chempathway$pathwayRampId)
  pathway_db <- pathway_db[!duplicated(tolower(names(pathway_db)))]

  if (methods::is(object, "SummarizedExperiment") && missing(slot)) {
    slot <- .integrationAssay(.experimentForAssay(object, assay))
  }
  E <- .assayData(object, assay, slot)

  pathway_sums <- list()
  for (i in seq_along(pathway_db)) {
    pathway <- pathway_db[[i]]
    pathway <- intersect(unique(pathway), rownames(E))
    if (!length(pathway)) {
      if (remove.nans) next
      score <- rep(NA_real_, ncol(object))
    } else {
      score <- .pathwayScores(list(pathway), object, assay, slot)[, 1L]
    }
    pathway_sums[[names(pathway_db)[i]]] <- score
  }

  if (!length(pathway_sums)) {
    stop("No pathways contain features from the selected assay.", call. = FALSE)
  }
  pathway_mtx <- do.call(cbind, pathway_sums)
  colnames(pathway_mtx) <- names(pathway_sums)

  pathway_mtx <- t(pathway_mtx)

  colnames(pathway_mtx) <- colnames(object)

  pathwayMetadata <- chempathway %>%
    filter(pathwayRampId %in% rownames(pathway_mtx)) %>%
    select(pathwayRampId, pathwayName) %>%
    distinct()
  rownames(pathwayMetadata) <- pathwayMetadata$pathwayRampId
  pathwayMetadata <- pathwayMetadata[rownames(pathway_mtx), , drop = FALSE]

  pathwayExperiment <- SingleCellExperiment::SingleCellExperiment(
    assays = list(pathwayScores = pathway_mtx),
    rowData = S4Vectors::DataFrame(pathwayMetadata)
  )
  SingleCellExperiment::altExp(object, new.assay) <- pathwayExperiment

  return(object)
}









############################# PATHWAY HELPER FUNCTIONS ########################################

#' Helper function for building a pathway db based on detected
#'
#' @param input_id Vector of characters defining the detected .
#' @param analytehaspathway A dataframe containing RAMP_pathway ID's.
#' @param chem_props A database containing the chemical properties and metadata of each RAMP_DB analyte.
#' @param pathway A dataframe containing RAMP_DB pathways and their relative metadata
#'
#' @return A analyte database containing corresponding pathways associated with each detected metabolite
#'
#' @examples
#' #HELPER FUNCTION
get_analytes_db <- function(input_id,analytehaspathway,chem_props,pathway) {

  rampid = unique(chem_props$ramp_id[which(chem_props$chem_source_id %in% unique(input_id))])
  #
  pathway_ids = unique(analytehaspathway$pathwayRampId[which(analytehaspathway$rampId %in% rampid)])

  analytes_db = lapply(pathway_ids, function(x) {
    content = analytehaspathway$rampId[which(analytehaspathway$pathwayRampId == x)]
    content = content[which(grepl(content, pattern = "RAMP_C"))]
    return(content)
  })
  analytes_db_name = unlist(lapply(pathway_ids, function(x) {
    name = pathway$pathwayName[which(pathway$pathwayRampId == x)]
    return(name)
  }))
  names(analytes_db) = analytes_db_name
  return(analytes_db)
}
