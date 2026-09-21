# Equivalence requires identical tested members in every supplied region.
# This is an evidence audit; it never changes the multiple-testing family.
.wf_pathway_equivalence <- function(regional) {
  coverage <- attr(regional, "pathway_coverage")
  ids <- unique(regional$pathwayRampId)
  if (!length(coverage) || !length(ids)) return(data.frame())
  signature <- vapply(ids, function(id) paste(vapply(coverage, function(c) {
    members <- c$used_members[[match(id, c$pathwayRampId)]]
    paste(sort(members), collapse = "|")
  }, character(1)), collapse = "::"), character(1))
  group <- match(signature, unique(signature))
  parts <- split(seq_along(ids), group)
  data.frame(member_group = paste0("set_", names(parts)), pathway_labels = lengths(parts),
    pathway_ids = I(lapply(parts, function(i) ids[i])),
    pathway_names = I(lapply(parts, function(i) regional$pathwayName[match(ids[i], regional$pathwayRampId)])),
    row.names = NULL)
}

.wf_regional_pathways <- function(result) {
  field <- result$region_analysis$field
  if (is.null(field)) return(result)
  ranking <- indices <- provenance <- list()
  for (name in names(result$analysis)) {
    index <- result$resources[[name]]$pathway_index
    cfg <- result$settings$modalities[[name]]$pathways
    if (is.null(index) || identical(cfg$regional, FALSE)) next
    assay <- result$analysis[[name]]$pathways$identity_assay
    if (is.null(assay)) {
      assay <- paste0("workflow_ids_", name)
      spec <- result$settings$modalities[[name]]
      result$object <- createPathwayAssay(result$object, analyte_type = "genes",
        assay = result$settings$assay_names[[name]], slot = "workflow", new_assay = assay,
        pathway_index = index, organism = spec$species, duplicate_genes = cfg$duplicate_genes %||% "error", verbose = FALSE)
      result$analysis[[name]]$pathways$identity_assay <- assay
    }
    markers <- findAllDEMs(result$object, ident = field, assay = assay, slot = "workflow",
      method = "markers", min_region_size = result$region_analysis$min_region_size, spatial_blocks = 0)
    if (!identical(markers$status, "completed")) next
    tab <- markers$DEMs
    mat <- vapply(colnames(markers$expression), function(region) {
      d <- tab[tab$cluster == region, ]
      d$logFC[match(rownames(markers$expression), d$gene)]
    }, numeric(nrow(markers$expression)))
    dimnames(mat) <- dimnames(markers$expression)
    type <- if (result$settings$modalities[[name]]$type == "metabolomics") "metabolites" else "genes"
    r <- stats::setNames(list(mat), type)
    enrichment <- findRegionalPathways(result$object, ident = field, ranks = r,
      pathway_index = index, min_path_size = cfg$min_size %||% 5, max_path_size = cfg$max_size %||% 500,
      verbose = FALSE)
    result$analysis[[name]]$pathways$regional <- enrichment
    result$analysis[[name]]$pathways$member_equivalence <- .wf_pathway_equivalence(enrichment)
    result$analysis[[name]]$pathways$identity_markers <- markers
    if (isTRUE(cfg$geseca) && any(result$analysis[[name]]$pathways$coverage$eligible)) {
      E <- .assayData(result$object, assay, "workflow")
      result$analysis[[name]]$pathways$geseca <- runRAMPGeseca(E, pathway_index = index,
        organism = result$settings$modalities[[name]]$species %||% "custom",
        minSize = cfg$min_size %||% 5, maxSize = cfg$max_size %||% 500,
        nproc = 1, nPermSimple = 1000, eps = 1e-10)
    }
    ranking[[name]] <- r; indices[[name]] <- index; provenance[[name]] <- index$provenance
  }
  # Combined enrichment is available only for a common index and disjoint
  # analyte types. Separate RNA assays or incompatible species are not merged.
  if (length(ranking) == 2L && length(unique(unlist(lapply(ranking, names)))) == 2L &&
      identical(provenance[[1]], provenance[[2]])) {
    r <- c(ranking[[1]], ranking[[2]])
    cfg <- result$settings$modalities[[names(ranking)[1]]]$pathways
    result$joint_pathways <- findRegionalPathways(result$object, ident = field, ranks = r,
      pathway_index = indices[[1]], min_path_size = cfg$min_size %||% 5, max_path_size = cfg$max_size %||% 500,
      verbose = FALSE)
  }
  result
}
