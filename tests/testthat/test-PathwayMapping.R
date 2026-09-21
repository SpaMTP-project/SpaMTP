pathwayClosureFixture <- function() {
  db <- list(source_df = data.frame(
    rampId = paste0("RAMP_G_", c(1, 2, 3, 4, 5, 5, 6)),
    sourceId = c("entrez:101", "uniprot:U1", "entrez:102", "entrez:103",
                 "entrez:101", "entrez:103", "entrez:104")),
    gene_reference = data.frame(hgnc_id = paste0("HGNC:", 1:4),
      symbol = c("GENEA", "GENEB", "GENEC", "GENED"), entrez_id = as.character(101:104),
      uniprot_ids = paste0("U", 1:4), prev_symbol = c("OLDA", "", "", "")),
    analytehaspathway = data.frame(rampId = paste0("RAMP_G_", c(1, 2, 3, 6, 2, 4, 5, 6)),
      pathwayRampId = c("P1", "P1", "P1", "P1", "P2", "P2", "conflict", "unmeasured")),
    pathway = data.frame(pathwayRampId = c("P1", "P2", "conflict", "unmeasured"),
      pathwayName = c("shared name", "shared name", "conflict", "unmeasured")))
  db$analyte <- data.frame(rampId = unique(db$source_df$rampId))
  object <- nativeFixture()
  E <- rbind(GENEA = c(1, 2, 4, 8, 16, 32), GENEB = c(1, 3, 2, 6, 5, 8),
             GENEC = c(3, 2, 1, 5, 7, 9))
  colnames(E) <- colnames(object)
  SingleCellExperiment::altExp(object, "transcriptome") <-
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = E, logcounts = log1p(E)))
  SingleCellExperiment::reducedDim(object, "PCA") <- cbind(1:6, c(1, 3, 2, 5, 4, 6))
  list(database = db, object = object, expression = E, index = buildPathwayIndex(db))
}

test_that("shared membership and coverage are identical in Fisher and expression scores", {
  f <- pathwayClosureFixture()
  idx <- f$index
  expect_setequal(idx$members$P1, c("RAMP_G_1", "RAMP_G_3", "RAMP_G_6"))
  expect_length(idx$raw_members$P1, 4L)
  expect_length(idx$members$conflict, 0L)
  expect_equal(idx$conflicting_members$conflict, "RAMP_G_5")
  original <- f$object
  mapped <- createPathwayAssay(original, analyte_type = "genes", assay = "transcriptome",
    slot = "logcounts", new_assay = "gene_ids", pathway_index = idx)
  E <- .assayData(mapped, "gene_ids", "logcounts")
  expect_identical(rownames(E), c("RAMP_G_1", "RAMP_G_3", "RAMP_G_4"))
  expect_equal(as.matrix(E), `rownames<-`(log1p(f$expression), rownames(E)))
  expect_identical(SummarizedExperiment::assayNames(SingleCellExperiment::altExp(mapped, "gene_ids")), "logcounts")
  direct <- createPathwayObject(original, assay = "transcriptome", slot = "logcounts", pathway_index = idx)
  indirect <- createPathwayObject(mapped, assay = "gene_ids", pathway_index = idx)
  scores <- .assayData(direct, "pathway", "pathwayScores")
  expect_equal(scores, .assayData(indirect, "pathway", "pathwayScores"))
  expected <- as.numeric(scale(colSums(log1p(f$expression[1:2, ])) / sqrt(2)))
  expect_equal(as.numeric(scores["P1", ]), expected)
  fisher <- fishersPathwayAnalysis(list(genes = "GENEA"),
    universe = list(genes = rownames(f$expression)), pathway_index = idx, min_path_size = 1,
    verbose = FALSE)
  coverage <- S4Vectors::metadata(SingleCellExperiment::altExp(direct, "pathway"))$pathway_mapping$coverage
  fcoverage <- attr(fisher, "pathway_coverage")
  expect_identical(coverage$measured_members, fcoverage$measured_members)
  expect_identical(coverage$database_members, fcoverage$database_members)
  p1 <- coverage[coverage$pathwayRampId == "P1", ]
  expect_equal(p1$database_size, 3L)
  expect_equal(p1$database_raw_size, 4L)
  expect_equal(p1$measured_size, 2L)
  expect_equal(p1$coverage_fraction, 2 / 3)
  expect_equal(p1$used_size, 2L)
  expect_equal(coverage$excluded_conflict_count[coverage$pathwayRampId == "conflict"], 1L)
  expect_true(is.na(coverage$coverage_fraction[coverage$pathwayRampId == "conflict"]))
  expect_equal(coverage$coverage_fraction[coverage$pathwayRampId == "unmeasured"], 0)
  expect_identical(SummarizedExperiment::assays(original), SummarizedExperiment::assays(direct))
  expect_identical(SingleCellExperiment::altExp(original, "transcriptome"), SingleCellExperiment::altExp(direct, "transcriptome"))
  expect_identical(SpatialExperiment::spatialCoords(original), SpatialExperiment::spatialCoords(direct))
})

test_that("symbol, previous symbol, stable ID and equivalent RaMP expression agree", {
  f <- pathwayClosureFixture()
  expected <- createPathwayObject(f$object, assay = "transcriptome", pathway_index = f$index)
  for (ids in list(c("OLDA", "GENEB", "GENEC"), c("101", "102", "103"),
                   c("RAMP_G_2", "RAMP_G_3", "RAMP_G_4"))) {
    x <- f$object
    rna <- SingleCellExperiment::altExp(x, "transcriptome")
    rownames(rna) <- ids
    SingleCellExperiment::altExp(x, "transcriptome") <- rna
    observed <- createPathwayObject(x, assay = "transcriptome", pathway_index = f$index)
    expect_equal(.assayData(observed, "pathway", "pathwayScores"), .assayData(expected, "pathway", "pathwayScores"))
  }
  x <- f$object
  rna <- SingleCellExperiment::altExp(x, "transcriptome")
  rna <- rna[c(1, 1, 2, 3), ]
  rownames(rna) <- c("GENEA", "OLDA", "GENEB", "GENEC")
  SingleCellExperiment::altExp(x, "transcriptome") <- rna
  duplicated <- createPathwayObject(x, assay = "transcriptome", pathway_index = f$index)
  expect_equal(.assayData(duplicated, "pathway", "pathwayScores"), .assayData(expected, "pathway", "pathwayScores"))
  SummarizedExperiment::assay(rna, "counts")[2, ] <- 3 * f$expression[1, ]
  SingleCellExperiment::altExp(x, "transcriptome") <- rna
  expect_error(createPathwayAssay(x, "genes", "transcriptome", pathway_index = f$index), "Multiple expression")
  averaged <- createPathwayAssay(x, "genes", "transcriptome", pathway_index = f$index, duplicate_genes = "mean")
  expect_equal(as.numeric(.assayData(averaged, "pathway")["RAMP_G_1", ]), 2 * as.numeric(f$expression[1, ]))
})

test_that("named plots use score-object members and reject ambiguous display names", {
  f <- pathwayClosureFixture()
  object <- createPathwayObject(f$object, assay = "transcriptome", pathway_index = f$index)
  expected <- .assayData(object, "pathway", "pathwayScores")
  reduced <- plotPathways(c("P1", "P2"), f$object, assay = "transcriptome", pathway_index = f$index)
  spatial <- plotPathwaysSpatially(c("P1", "P2"), f$object, assay = "transcriptome", pathway_index = f$index)
  for (id in c("P1", "P2")) expect_identical(spatial[[id]]$data$score, unname(expected[id, ]))
  for (plots in list(reduced, spatial)) {
    expect_identical(attr(plots, "pathway_coverage")$used_size, c(2L, 2L))
    for (id in c("P1", "P2")) {
      expect_equal(attr(plots[[id]], "pathway_scores"), expected[id, ])
      expect_s3_class(plots[[id]], "ggplot")
      expect_true(length(ggplot2::ggplot_build(plots[[id]])$data) > 0)
    }
  }
  expect_error(plotPathways("shared name", f$object, assay = "transcriptome", pathway_index = f$index), "ambiguous")
})

test_that("coverage records zero overlap and size exclusions without changing the universe", {
  f <- pathwayClosureFixture()
  result <- fishersPathwayAnalysis(list(genes = "GENEC"), universe = list(genes = rownames(f$expression)),
    pathway_index = f$index, min_path_size = 1, verbose = FALSE)
  expect_equal(result$p_val[result$pathwayRampId == "P1"], 1)
  expect_equal(result$background_analytes_number, c(3L, 3L))
  scored <- createPathwayObject(f$object, assay = "transcriptome", pathway_index = f$index, remove.nans = FALSE)
  expect_true(all(is.na(.assayData(scored, "pathway", "pathwayScores")["unmeasured", ])))
  expect_error(createPathwayObject(f$object, assay = "transcriptome", pathway_index = f$index,
                                   min_path_size = 3), "No pathways")
  changed <- f$database
  changed$analytehaspathway$pathwayRampId[1] <- "P2"
  expect_error(createPathwayObject(f$object, assay = "transcriptome", database = changed,
                                   pathway_index = f$index), "does not match")
  x <- f$object
  S4Vectors::metadata(x)$SpaMTPData <- list(organism = "Mus musculus")
  expect_error(createPathwayObject(x, assay = "transcriptome", pathway_index = f$index), "species")
})

test_that("GESECA receives canonical expression and pathway IDs instead of display names", {
  f <- pathwayClosureFixture()
  f$database$analytehaspathway <- rbind(f$database$analytehaspathway,
    data.frame(rampId = "RAMP_C_1", pathwayRampId = "P1"))
  f$index <- buildPathwayIndex(f$database)
  seen <- NULL
  local_mocked_bindings(geseca = function(pathways, E, ...) {
    seen <<- list(pathways = pathways, E = E)
    data.frame(pathway = names(pathways))
  }, .package = "fgsea")
  result <- runRAMPGeseca(f$expression, pathway_index = f$index)
  expect_setequal(rownames(seen$E), c("RAMP_G_1", "RAMP_G_3", "RAMP_G_4"))
  expect_setequal(names(seen$pathways), c("P1", "P2"))
  expect_setequal(seen$pathways$P1, c("RAMP_G_1", "RAMP_G_3"))
  expect_true("coverage_fraction" %in% names(attr(result, "pathway_coverage")))
  coverage <- attr(result, "pathway_coverage")
  expect_equal(coverage$database_size[coverage$pathwayRampId == "P1"], 3L)
})

test_that("compound-only indices without crossreferences are reusable in Fisher", {
  db <- list(analytehaspathway = data.frame(rampId = c("RAMP_C_1", "RAMP_C_2"),
    pathwayRampId = c("P1", "P2")),
    pathway = data.frame(pathwayRampId = c("P1", "P2", "empty")))
  index <- buildPathwayIndex(db)
  result <- fishersPathwayAnalysis(list(metabolites = "RAMP_C_1"),
    universe = list(metabolites = c("RAMP_C_1", "RAMP_C_2")), pathway_index = index,
    min_path_size = 1, verbose = FALSE)
  expect_equal(result$p_val, c(.5, 1))
  coverage <- attr(result, "pathway_coverage")
  expect_equal(coverage$used_size, c(1L, 1L, 0L))
})

test_that("RaMP-only assays exclude unmapped symbols from the measured universe", {
  f <- pathwayClosureFixture()
  db <- list(source_df = data.frame(rampId = c("RAMP_G_1", "RAMP_G_2"),
    sourceId = c("entrez:101", "entrez:102"), commonName = c("GENEA", "GENEB")),
    analytehaspathway = data.frame(rampId = c("RAMP_G_1", "RAMP_G_2"),
      pathwayRampId = c("P1", "P2")), pathway = data.frame(pathwayRampId = c("P1", "P2")))
  index <- buildPathwayIndex(db)
  expect_warning(mapped <- createPathwayAssay(f$object, "genes", "transcriptome",
    slot = "logcounts", pathway_index = index), "unresolved in RaMP-only")
  target <- SingleCellExperiment::altExp(mapped, "pathway")
  expect_identical(rownames(target), c("RAMP_G_1", "RAMP_G_2"))
  expect_equal(S4Vectors::metadata(target)$pathway_mapping$inputs$status,
    c("mapped_ramp_only", "mapped_ramp_only", "unmapped"))
  fisher <- suppressWarnings(fishersPathwayAnalysis(list(genes = "GENEA"),
    universe = list(genes = rownames(f$expression)), pathway_index = index,
    min_path_size = 1, verbose = FALSE))
  expect_identical(S4Vectors::metadata(target)$pathway_mapping$coverage$measured_members,
    attr(fisher, "pathway_coverage")$measured_members)
})
