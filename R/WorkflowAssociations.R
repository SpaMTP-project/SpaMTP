.wf_associations <- function(result) {
  maximum <- result$settings$association_features
  names <- names(result$analysis)
  if (length(names) < 2L || !maximum) return(result)
  fields <- unique(c("sample_id", result$settings$replicate, result$region_analysis$field))
  md <- as.data.frame(SummarizedExperiment::colData(result$object))
  fields <- fields[fields %in% colnames(md)]
  keep <- stats::complete.cases(md[, fields, drop = FALSE])
  for (f in fields) if (is.character(md[[f]]) || is.factor(md[[f]]))
    keep <- keep & nzchar(trimws(as.character(md[[f]])))
  object <- result$object[, keep, drop = FALSE]
  if (ncol(object) < 4L) return(result)
  answers <- list()
  for (pair in utils::combn(names, 2, simplify = FALSE)) {
    a <- pair[1]; b <- pair[2]
    # Anchors come from the regional marker ranking, round-robin by region.
    # Without regions, the explicitly recorded variable-feature screen is used.
    markers <- result$region_analysis$modalities[[a]]$tables
    ranked <- lapply(markers, function(t) t$feature[order(-t$mean_auc)])
    anchors <- unique(unlist(lapply(seq_len(maximum), function(i)
      vapply(ranked, function(x) if (length(x) >= i) x[i] else NA_character_, character(1)))))
    anchors <- head(unique(c(stats::na.omit(anchors), result$analysis[[a]]$pca$features)), maximum)
    targets <- result$analysis[[b]]$pca$features
    if (!length(anchors) || !length(targets)) next
    # The public interface accepts generic feature IDs in its historical mz/gene
    # arguments. A temporary paired view allows any two of three or more omics.
    ea <- .experimentForAssay(object, result$settings$assay_names[[a]])
    paired <- SpatialExperiment::SpatialExperiment(
      assays = list(workflow = SummarizedExperiment::assay(ea, "workflow")),
      rowData = SummarizedExperiment::rowData(ea), colData = SummarizedExperiment::colData(object),
      spatialCoords = SpatialExperiment::spatialCoords(object))
    SingleCellExperiment::altExp(paired, "target") <-
      .experimentForAssay(object, result$settings$assay_names[[b]])
    query <- lapply(anchors, function(anchor) {
      t <- findCorrelatedFeatures(paired, mz = anchor, ST.assay = "target",
        SM.slot = "workflow", ST.slot = "workflow", nfeatures = NULL,
        covariates = fields, features = list(metabolite = anchor, gene = targets))
      t <- t[t$modality == "gene", , drop = FALSE]
      data.frame(feature1 = anchor, feature2 = t$features,
        correlation = t$correlation, correlation_raw = t$correlation_raw,
        modality1 = a, modality2 = b, observations = ncol(paired))
    })
    t <- do.call(rbind, query)
    t <- t[order(-abs(t$correlation), t$feature1, t$feature2, na.last = TRUE), , drop = FALSE]
    attr(t, "association") <- list(function_name = "findCorrelatedFeatures",
      covariates = fields, anchors = anchors, targets = targets,
      observations = colnames(object), excluded_observations = colnames(result$object)[!keep],
      selection = "Anchors: regional AUC rank, round-robin across regions, then variable features; targets: all PCA-selected features",
      interpretation = "Raw and covariate-residual Pearson correlations on the same observations. Descriptive screening; no pixel-based P values. Selection uses these data and is not independent validation.")
    answers[[paste(pair, collapse = "__")]] <- t
  }
  result$associations <- answers
  result$stages$detail[result$stages$stage == "associations"] <-
    "findCorrelatedFeatures: regional marker anchors; raw versus sample/replicate/region-residual correlations on matched observations."
  result$stages$status[result$stages$stage == "associations"] <- if (length(answers)) "completed" else "skipped"
  result
}
