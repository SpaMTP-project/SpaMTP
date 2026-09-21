singleFeatureMarkerFixture <- function(sparse = FALSE) {
  withr::local_seed(1921)
  group <- rep(rep(c("A", "B"), each = 3), 4)
  values <- log1p(stats::rpois(24, ifelse(group == "A", 6, 16)))
  E <- matrix(values, nrow = 1,
    dimnames = list("feature1", paste0("p", seq_along(values))))
  if (sparse) E <- Matrix::Matrix(E, sparse = TRUE)
  SpatialExperiment::SpatialExperiment(assays = list(logcounts = E),
    colData = S4Vectors::DataFrame(region = group,
      patient = rep(paste0("patient", 1:4), each = 6)),
    spatialCoords = cbind(x = rep(1:6, 4), y = rep(1:4, each = 6)))
}

test_that("single-feature markers preserve dimensions with spatial sensitivity", {
  for (sparse in c(FALSE, TRUE)) for (blocks in c(0, 4)) {
    object <- singleFeatureMarkerFixture(sparse)
    before <- object
    fit <- findAllDEMs(object, "region", assay = "main", slot = "logcounts",
      method = "markers", spatial_blocks = blocks)
    E <- as.matrix(SummarizedExperiment::assay(object))
    expected <- matrix(c(mean(E[, object$region == "A"]),
      mean(E[, object$region == "B"])), 1,
      dimnames = list("feature1", c("A", "B")))
    expect_identical(fit$status, "completed")
    expect_equal(fit$expression, expected)
    expect_identical(dim(fit$expression), c(1L, 2L))
    expect_equal(fit$DEMs$gene, rep("feature1", 2))
    expect_true(all(is.finite(fit$DEMs$mean_auc)))
    expect_equal(fit$DEMs$detected_reference, rep(1, 2))
    if (blocks) {
      expect_true(all(fit$DEMs$stability_fits > 0))
      expect_equal(fit$DEMs$rank_best_block, rep(1, 2))
      expect_equal(fit$DEMs$rank_worst_block, rep(1, 2))
      expect_true(all(is.finite(fit$DEMs$direction_stability)))
    }
    expect_identical(object, before)
  }
})

test_that("single-feature replicate summaries and contrasts remain matrices", {
  for (sparse in c(FALSE, TRUE)) for (paired in c(FALSE, TRUE)) {
    object <- singleFeatureMarkerFixture(sparse)
    if (!paired) object$patient <- paste(object$patient, object$region, sep = "_")
    summary <- findAllDEMs(object, "region", assay = "main", slot = "logcounts",
      method = "replicate")
    expect_identical(dim(summary$group_means), c(1L, 2L))
    expect_identical(rownames(summary$group_means), "feature1")
    fit <- findAllDEMs(object, "region", assay = "main", slot = "logcounts",
      method = "replicate", replicate = "patient",
      contrasts = list(B_A = list(numerator = "B", denominator = "A", paired = paired)))
    result <- fit$tests$B_A
    expect_identical(result$status, "completed")
    expect_identical(dim(result$expression), c(1L, 8L))
    expect_identical(rownames(result$expression), "feature1")
    E <- as.matrix(SummarizedExperiment::assay(object))
    expected <- matrix(vapply(seq_len(nrow(result$units)), function(i) mean(
      E[, object$patient == result$units$replicate[i] &
        object$region == result$units$group[i]]), numeric(1)), 1,
      dimnames = dimnames(result$expression))
    expect_equal(result$expression, expected)
    reference <- limma::lmFit(expected, result$design)
    reference <- limma::eBayes(limma::contrasts.fit(reference, result$coefficient),
      trend = TRUE)
    table <- limma::topTable(reference, number = Inf, sort.by = "none", confint = TRUE)
    expect_equal(result$table$logFC, table$logFC)
    expect_equal(result$table$P.Value, table$P.Value)
    expect_true(all(is.finite(result$table$P.Value)))
  }
})
