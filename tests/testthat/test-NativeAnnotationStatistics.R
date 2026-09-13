annotationStatisticsFixture <- function() {
  object <- nativeFixture()
  SummarizedExperiment::rowData(object)$all_Ramp_IDs <-
    c("RAMP_C_1; RAMP_C_2; RAMP_C_3; RAMP_C_1", "RAMP_C_1", "RAMP_C_1; RAMP_C_2")
  pathways <- rbind(pA = 1:6, pB = c(2, 1, 5, 3, 6, 4), pC = 6:1, pD = rep(2, 6))
  colnames(pathways) <- colnames(object)
  SingleCellExperiment::altExp(object, "scores") <-
    SingleCellExperiment::SingleCellExperiment(
      assays = list(pathwayScores = pathways, negative = -pathways))
  database <- list(
    analytehaspathway = data.frame(
      rampId = c("RAMP_C_1", "RAMP_C_1", "RAMP_C_2", "RAMP_C_3"),
      pathwayRampId = c("pA", "pB", "pB", "pC")),
    source_df = data.frame(rampId = paste0("RAMP_C_", 1:3),
      sourceId = paste0("id", 1:3), commonName = c("First", "Second", "Third")),
    pathway = data.frame(pathwayRampId = c("pA", "pB", "pC"),
      pathwayName = c("A", "B", "C")))
  list(object = object, database = database)
}

test_that("native annotation ranking preserves the finite legacy formula", {
  fixture <- annotationStatisticsFixture()
  object <- fixture$object
  before <- object
  result <- calculateSingleAnnotationStatistics(100, object,
    pathway.assay = "scores", database = fixture$database)
  correlation <- stats::cor(1:6, c(2, 1, 5, 3, 6, 4))
  maximum <- c(1, correlation, -1)
  number <- c(2, 1, 0)
  expected <- as.numeric(scale(abs(maximum))) + as.numeric(scale(number))
  result <- result[match(paste0("RAMP_C_", 1:3), result$ramp_id), ]
  expect_equal(result$max_cor, maximum)
  expect_equal(result$n_sig_path, number)
  expect_equal(result$z_score, expected)
  expect_equal(result$pval, stats::pnorm(expected, lower.tail = FALSE))
  expect_equal(result$pval_adj, stats::p.adjust(result$pval, "BH"))
  expect_identical(result$metabolite, c("First", "Second", "Third"))
  reordered <- calculateSingleAnnotationStatistics("a", object[, 6:1],
    pathway.assay = "scores", database = fixture$database)
  expect_equal(reordered$z_score[match(result$ramp_id, reordered$ramp_id)], result$z_score)
  negative <- calculateSingleAnnotationStatistics("a", object,
    pathway.assay = "scores", pathway.slot = "negative", database = fixture$database)
  expect_equal(negative$max_cor[match(result$ramp_id, negative$ramp_id)], -maximum)
  expect_identical(object, before)
  expect_error(calculateSingleAnnotationStatistics("a", object, pathway.assay = "scores",
    pathway.slot = "absent", database = fixture$database), "was not found")
})

test_that("constant, absent and duplicated candidates do not create invalid ranks", {
  fixture <- annotationStatisticsFixture()
  object <- fixture$object
  result <- calculateSingleAnnotationStatistics("constant", object,
    pathway.assay = "scores", database = fixture$database)
  expect_true(all(is.na(result$max_cor)))
  expect_true(all(is.na(result$z_score)))
  expect_true(all(is.na(result$pval_adj)))
  expect_error(calculateSingleAnnotationStatistics("b", object,
    pathway.assay = "scores", database = fixture$database), "two candidate")
  expect_error(calculateSingleAnnotationStatistics("a", object, corr_weight = -1,
    pathway.assay = "scores", database = fixture$database), "Weights")
  expect_error(calculateSingleAnnotationStatistics("a", object, corr_theshold = 2,
    pathway.assay = "scores", database = fixture$database), "between")
  fixture$database$analytehaspathway$pathwayRampId[] <- "pA"
  equal <- calculateSingleAnnotationStatistics("a", object,
    pathway.assay = "scores", database = fixture$database)
  expect_equal(equal$z_score, c(0, 0, 0))
  expect_equal(equal$pval_adj, c(0.5, 0.5, 0.5))
})

test_that("batch ranking reuses pathway scores and preserves feature metadata", {
  fixture <- annotationStatisticsFixture()
  object <- fixture$object
  top <- calculateAnnotationStatistics(object, pathway.assay = "scores",
    pathway.slot = "pathwayScores", pathway.scores = TRUE, database = fixture$database)
  expect_identical(top$mz_names, rownames(object))
  expect_identical(top$annotation, letters[1:3])
  expect_identical(top$mz, c(100, 200, 300))
  expect_identical(top$ramp_id[1], "RAMP_C_1")
  expect_true(all(is.na(top$ramp_id[2:3])))
  all <- calculateAnnotationStatistics(object, pathway.assay = "scores",
    pathway.slot = "pathwayScores", pathway.scores = TRUE, return.top = FALSE,
    database = fixture$database)
  expect_identical(names(all), rownames(object))
  expect_null(all$b)
  analytes <- rbind(RAMP_C_1 = 1:6, RAMP_C_2 = c(2, 1, 5, 3, 6, 4), RAMP_C_3 = 6:1)
  colnames(analytes) <- colnames(object)
  SingleCellExperiment::altExp(object, "analytes") <-
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = analytes))
  scored <- createPathwayObject(object, assay = "analytes", slot = "counts",
    new.assay = "expected", database = fixture$database)
  expected <- calculateAnnotationStatistics(scored, pathway.assay = "expected",
    pathway.slot = "pathwayScores", pathway.scores = TRUE, database = fixture$database)
  actual <- calculateAnnotationStatistics(object, pathway.assay = "analytes",
    database = fixture$database)
  expect_equal(actual, expected)
  expect_false(".annotationPathways" %in% SingleCellExperiment::altExpNames(object))
})

test_that("current RaMP candidates take precedence over legacy rowData", {
  fixture <- annotationStatisticsFixture()
  object <- fixture$object
  current <- data.frame(observed_mz = 100, Adduct = "M+H",
    Ramp_IDs = "RAMP_C_2; RAMP_C_3", Score = 1, MassScore = 1,
    ChemicalScore = 1, IsotopeScore = NA_real_, AdductNetworkScore = NA_real_)
  object <- .setStoredData(object, "mz_annotation",
    list(results = current, metadata = list(schema_version = 2L)))
  result <- calculateSingleAnnotationStatistics("a", object,
    pathway.assay = "scores", database = fixture$database)
  expect_setequal(result$ramp_id, c("RAMP_C_2", "RAMP_C_3"))
  object <- fixture$object
  SummarizedExperiment::rowData(object)$all_Ramp_IDs <- NULL
  SummarizedExperiment::rowData(object)$all_Isomers_IDs <- c("id1; id2; id3", "id1", "id1; id2")
  result <- calculateSingleAnnotationStatistics("a", object,
    pathway.assay = "scores", database = fixture$database)
  expect_setequal(result$ramp_id, paste0("RAMP_C_", 1:3))
})

test_that("Pearson correlation agrees with Cardinal colocalization", {
  object <- nativeFixture()
  cardinal <- asCardinal(object, assayName = "counts")
  result <- suppressWarnings(Cardinal::colocalized(cardinal, mz = 100,
    verbose = FALSE, BPPARAM = BiocParallel::SerialParam()))
  observed <- result$cor[match(200, result$mz)]
  expect_equal(observed, stats::cor(1:6, c(2, 1, 5, 3, 6, 4)), tolerance = 1e-6)
})
