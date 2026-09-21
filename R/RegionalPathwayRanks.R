.regionalPathwaysFromRanks <- function(ranks, index, min_size, max_size, rank_scale, nPermSimple) {
  .pathway_validate_index(index)
  .validateFeatureCount(nPermSimple, "nPermSimple")
  if (!is.list(ranks) || !length(ranks) || is.null(names(ranks)) ||
      anyDuplicated(names(ranks)) || any(!names(ranks) %in% c("genes", "metabolites")))
    stop("ranks must be a named list of genes/metabolites feature-by-region matrices.", call. = FALSE)
  mapped <- list(); mapping <- list(); regions <- NULL
  for (type in names(ranks)) {
    E <- ranks[[type]]
    if (length(dim(E)) != 2L || !is.numeric(E) || is.null(rownames(E)) ||
        is.null(colnames(E)) || anyNA(rownames(E)) || anyNA(colnames(E)) ||
        anyDuplicated(rownames(E)) || anyDuplicated(colnames(E)))
      stop("Rank matrices require distinct feature and region names.", call. = FALSE)
    if (is.null(regions)) regions <- colnames(E)
    if (!setequal(regions, colnames(E))) stop("Rank matrices must describe the same regions.", call. = FALSE)
    if (type == "metabolites" && !all(grepl("^RAMP_C_", rownames(E))))
      stop("Metabolite ranks require explicit RaMP compound identities; use createPathwayAssay first.", call. = FALSE)
    m <- .pathway_expression(E[, regions, drop = FALSE], index, "error")
    ids <- rownames(m$expression)
    prefix <- if (type == "genes") "^RAMP_G_" else "^RAMP_C_"
    if (any(!grepl(prefix, ids))) stop("Rank identities do not match analyte type.", call. = FALSE)
    mapped[[type]] <- as.matrix(m$expression)
    mapping[[type]] <- m$mapping
  }
  measured <- unique(unlist(lapply(mapped, rownames)))
  coverage <- used_ranks <- tables <- list()
  for (region in regions) {
    blocks <- lapply(mapped, function(E) {
      v <- stats::setNames(E[, region], rownames(E))
      v <- v[is.finite(v)]
      if (rank_scale == "rms" && length(v)) {
        divisor <- sqrt(mean(v^2)); if (divisor > 0) v <- v / divisor
      }
      v
    })
    stat <- unlist(blocks, use.names = FALSE)
    names(stat) <- unlist(lapply(blocks, names), use.names = FALSE)
    if (anyDuplicated(names(stat))) stop("An identity was repeated across rank blocks.", call. = FALSE)
    used_ranks[[region]] <- sort(stat, decreasing = TRUE)
    c <- .pathway_coverage(index, measured, used_ids = names(stat),
      types = ifelse(names(ranks) == "genes", "G", "C"), min_size = min_size, max_size = max_size)
    coverage[[region]] <- c
    sets <- .pathway_sets_from_coverage(c)
    if (!length(sets) || length(stat) < 2L || all(stat == 0)) next
    t <- as.data.frame(fgsea::fgseaMultilevel(pathways = sets, stats = used_ranks[[region]],
      minSize = min_size, maxSize = max_size, nPermSimple = nPermSimple, eps = 1e-10, nproc = 1))
    if (!nrow(t)) next
    # Include all eligible pathways in the declared adjustment family, including
    # pathways for which the Monte Carlo algorithm could not estimate a P value.
    t$padj <- stats::p.adjust(t$pval, method = "BH", n = length(sets))
    names(t)[names(t) == "pathway"] <- "pathwayRampId"
    t$Cluster_id <- region
    t$tested_family <- length(sets)
    t <- merge(t, index$metadata, by = "pathwayRampId", all.x = TRUE, sort = FALSE)
    tables[[region]] <- t
  }
  answer <- if (length(tables)) do.call(rbind, tables) else data.frame(
    pathwayRampId = character(), pathwayName = character(), Cluster_id = character(),
    NES = numeric(), pval = numeric(), padj = numeric(), size = integer())
  rownames(answer) <- NULL
  attr(answer, "pathway_coverage") <- coverage
  attr(answer, "pathway_index") <- index$provenance
  attr(answer, "gene_mapping") <- mapping
  attr(answer, "ranks") <- used_ranks
  attr(answer, "rank_analysis") <- list(input = "complete unfiltered feature rankings",
    rank_scale = rank_scale, analyte_types = names(ranks), nPermSimple = nPermSimple,
    family = "All measured/used-size-eligible pathways, BH within region and configured modality combination",
    interpretation = "Competitive feature-set enrichment conditional on supplied ranks; P values do not test between-specimen regional reproducibility. Joint ranks are RMS-scaled per modality; pathways with many measured genes can still be gene-dominated.")
  answer
}
