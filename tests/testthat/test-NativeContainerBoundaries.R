test_that("normalization targets the selected experiment and validates intensities", {
  object <- nativeFixture()
  original <- SummarizedExperiment::assays(object)
  result <- normalizeSMData(object, "LogNormalize", assay = "transcriptome",
                            scale.factor = 100, verbose = FALSE)
  expect_equal(SummarizedExperiment::assays(result), original)
  transcriptome <- SingleCellExperiment::altExp(result, "transcriptome")
  counts <- SummarizedExperiment::assay(transcriptome, "counts")
  expected <- log1p(sweep(counts, 2, colSums(counts), "/") * 100)
  expect_equal(as.matrix(SummarizedExperiment::assay(transcriptome, "logcounts")),
               expected)
  expect_error(normalizeSMData(object, scale.factor = 0), "scale.factor")
  SummarizedExperiment::assay(object, "counts")[1, 1] <- -1
  expect_error(normalizeSMData(object), "non-negative")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(c(0, 0, 1, 3), 2)))
  sce <- normalizeSMData(sce, scale.factor = 10, verbose = FALSE)
  expect_equal(as.numeric(colSums(SummarizedExperiment::assay(sce, "normcounts"))),
               c(0, 10))
})

test_that("TMM supports a plain alternative SummarizedExperiment", {
  object <- nativeFixture()
  transcriptome <- SingleCellExperiment::altExp(object, "transcriptome")
  plain <- SummarizedExperiment::SummarizedExperiment(
    assays = SummarizedExperiment::assays(transcriptome))
  object <- addTranscriptome(object, plain)
  result <- tmmNormalize(object, "region", "edge", assay = "transcriptome")
  expect_equal(SummarizedExperiment::assays(result), SummarizedExperiment::assays(object))
  expect_true("normcounts" %in% SummarizedExperiment::assayNames(
    SingleCellExperiment::altExp(result, "transcriptome")))
  expect_error(tmmNormalize(object, "region", "edge", normalisation.type = "bad"),
                "arg")
})

test_that("paired pixels must be unique and SCE conversion retains infrastructure", {
  object <- nativeFixture()
  transcriptome <- .assayData(object, "transcriptome")
  colnames(transcriptome)[2] <- colnames(transcriptome)[1]
  expect_error(addTranscriptome(object, transcriptome), "unique")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = SummarizedExperiment::assays(object),
    rowData = SummarizedExperiment::rowData(object),
    colData = SummarizedExperiment::colData(object))
  sce$x <- 1:6
  sce$y <- rep(0, 6)
  SingleCellExperiment::altExps(sce) <- SingleCellExperiment::altExps(object)
  SingleCellExperiment::reducedDim(sce, "PCA") <- cbind(1:6, 6:1)
  SingleCellExperiment::colPair(sce, "neighbors") <-
    S4Vectors::SelfHits(1:5, 2:6, nnode = 6)
  result <- asSpatialExperiment(sce)
  expect_equal(SingleCellExperiment::altExps(result), SingleCellExperiment::altExps(sce))
  expect_equal(SingleCellExperiment::reducedDims(result), SingleCellExperiment::reducedDims(sce))
  expect_equal(SingleCellExperiment::colPairs(result), SingleCellExperiment::colPairs(sce))
  sce$horizontal <- sce$x
  sce$vertical <- sce$y
  renamed <- asSpatialExperiment(sce, coordinates = c("horizontal", "vertical"))
  expect_identical(colnames(SpatialExperiment::spatialCoords(renamed)), c("x", "y"))
  expect_s3_class(plotSpatialFeature(renamed, "a"), "ggplot")
  expect_error(asSpatialExperiment(sce, coordinates = "x"), "two distinct")
})

test_that("aligned binning supports mz-only rowData and conserves summed intensity", {
  object <- nativeFixture()
  SummarizedExperiment::rowData(object)$mz <- c(100, 100.5, 101.75)
  for (units in c("mz", "ppm")) {
    result <- binSpaMTP(object, resolution = if (units == "mz") 1 else 10000,
                         units = units)
    expect_equal(as.numeric(colSums(SummarizedExperiment::assay(result))),
                   as.numeric(colSums(SummarizedExperiment::assay(object))))
    expect_true(all(diff(SummarizedExperiment::rowData(result)$mz) > 0))
  }
  result <- binSpaMTP(object, resolution = 1, units = "mz")
  expect_equal(SummarizedExperiment::rowData(result)$mz, c(100, 101, 102))
  expect_error(binSpaMTP(object, resolution = 0), "resolution")
  expect_error(binSpaMTP(object, resolution = 10, assay = "transcriptome"), "primary")
  for (method in c("sum", "mean", "min", "max")) {
    result <- binSpaMTP(object, resolution = 10, units = "mz", method = method)
    counts <- SummarizedExperiment::assay(object)
    expected <- apply(counts, 2, switch(method, sum = sum, mean = mean, min = min, max = max))
    expect_equal(as.numeric(SummarizedExperiment::assay(result)), as.numeric(expected))
  }
})

test_that("integration excludes derived assays and standardizes merged modalities", {
  object <- nativeFixture()
  object <- createMergedModalityAssay(object, c("main", "transcriptome"))
  merged <- as.matrix(.assayData(object, "merged", "scaled"))
  expect_equal(as.numeric(rowMeans(merged)), rep(0, nrow(merged)), tolerance = 1e-12)
  expect_true(all(is.finite(merged)))
  result <- suppressWarnings(multiOmicIntegration(object,
    dims.list = list(1, 1), return.intermediate = TRUE))
  expect_identical(S4Vectors::metadata(result)$spamtp_integration$modalities,
                     c("main", "transcriptome"))
  expect_equal(ncol(SingleCellExperiment::reducedDim(result, "integrated")), 2)
  expect_true("merged" %in% SingleCellExperiment::altExpNames(result))
  expect_error(multiOmicIntegration(object, modalities = c("main", "absent")),
                "not found")
})

test_that("spatial graphs and plots separate interleaved tissue samples", {
  object <- nativeFixture(rep(c("s1", "s2"), 3))
  SpatialExperiment::spatialCoords(object) <-
    cbind(x = rep(1:3, each = 2), y = 0)
  result <- runSpatialGraphPCA(object, n_components = 2, slot = "counts",
                                n_neighbors = 2, platform = "ST", verbose = FALSE)
  graph <- S4Vectors::metadata(result)$spatialGraphs$SpatialKNN
  expect_true(all(graph[object$sample_id == "s1", object$sample_id == "s2"] == 0))
  expect_true(all(rowSums(graph) == 2))
  expect_equal(SingleCellExperiment::colPair(result, "SpatialKNN") |> length(), 12L)
  subset <- result[, 1:3]
  expect_equal(S4Vectors::nnode(SingleCellExperiment::colPair(subset, "SpatialKNN")), 3L)
  plot <- plotSpatialFeature(object, "a")
  points <- ggplot2::ggplot_build(plot)$data[[1]]
  expect_equal(length(unique(points$PANEL)), 2L)
})

test_that("Cardinal normalization never silently ignores an argument", {
  object <- asCardinal(nativeFixture())
  expect_error(normalizeSMData(object, "RC"), "Cardinal input supports TIC")
  expect_error(normalizeSMData(object, scale.factor = 100), "scale.factor")
  result <- binSpaMTP(object, resolution = 100, units = "mz")
  expect_s4_class(result, "SpatialExperiment")
  reference <- Cardinal::bin(object, resolution = 100, units = "mz", method = "sum")
  expect_equal(as.matrix(SummarizedExperiment::assay(result)),
                 as.matrix(Cardinal::spectra(reference)), ignore_attr = TRUE)
})

test_that("ROI app and standard subsetting use native metadata", {
  object <- nativeFixture()
  expect_s3_class(selectROIs(object, launch = FALSE), "shiny.appobj")
  result <- subsetSPM(object, subset = region == "edge", features = c("a", "b"))
  expect_identical(colnames(result), colnames(object)[c(1, 3, 5)])
  expect_identical(rownames(result), c("a", "b"))
  expect_identical(colnames(SingleCellExperiment::altExp(result, "transcriptome")),
                     colnames(result))
  expect_error(subsetSPM(object, cells = "absent"), "existing pixels")
  multi <- nativeFixture(rep(c("s1", "s2"), 3))
  expect_error(selectROIs(multi, launch = FALSE), "sampleId")
  expect_s3_class(selectROIs(multi, sampleId = "s2", launch = FALSE), "shiny.appobj")
})

test_that("technical pooling uses edgeR and limma treatment tables are valid", {
  withr::local_seed(29)
  counts <- matrix(stats::rpois(80 * 24, lambda = 20), 80, 24,
    dimnames = list(paste0("mz", 1:80), paste0("p", 1:24)))
  counts[1:5, 1:12] <- counts[1:5, 1:12] + 100
  object <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(group = rep(c("A.+_one", "B_two"), each = 12)))
  before <- .Random.seed
  pooled <- runPooling(object, "group", n = 3, assay = "main", slot = "counts",
                         verbose = FALSE)
  expect_identical(.Random.seed, before)
  expect_equal(ncol(pooled), 6L)
  expect_equal(as.numeric(rowSums(SummarizedExperiment::assay(pooled))),
                 as.numeric(rowSums(counts)))
  expect_equal(as.numeric(table(pooled$group)), c(3, 3))
  result <- runDE(pooled, object, "group", output_dir = NULL, run_name = "test",
                   n = 3, logFC_threshold = 1.2, annotation.column = NULL,
                   assay = "main", verbose = FALSE)
  expect_true(all(c("FDR", "regulate", "cluster") %in% colnames(result$DEMs)))
  expect_true(all(result$DEMs$FDR >= 0 & result$DEMs$FDR <= 1))
  strong <- result$DEMs[result$DEMs$gene == "mz1" & result$DEMs$cluster == "A.+_one", ]
  expect_equal(strong$regulate, "Up")
})

test_that("RaMP aggregation retains feature IDs and pixel ordering", {
  object <- nativeFixture()
  matrix <- rbind(geneA = 1:6, geneB = 6:1, geneC = c(3, 1, 4, 1, 5, 9))
  colnames(matrix) <- colnames(object)
  object <- addTranscriptome(object, matrix)
  database <- list(chem_props = data.frame(), source_df = data.frame(
    rampId = c("RAMP_G_2", "RAMP_G_1", "RAMP_G_2"),
    commonName = c("GENEA", "GENEC", "GENEB")))
  result <- createPathwayAssay(object, analyte_type = "genes", assay = "transcriptome",
                                new_assay = "ramp", database = database, verbose = FALSE)
  ramp <- .assayData(result, "ramp")
  expect_setequal(rownames(ramp), c("RAMP_G_1", "RAMP_G_2"))
  expect_identical(colnames(ramp), colnames(object))
  expect_equal(as.numeric(ramp["RAMP_G_2", ]), rep(3.5, 6))
  expect_equal(as.numeric(ramp["RAMP_G_1", ]), as.numeric(matrix["geneC", ]))
  pathwayDb <- list(analytehaspathway = data.frame(
    pathwayRampId = c("P1", "P2"), rampId = c("RAMP_G_2", "absent")),
    pathway = data.frame(pathwayRampId = c("P1", "P2"), pathwayName = c("one", "two")))
  result <- createPathwayObject(result, assay = "ramp", database = pathwayDb)
  scores <- .assayData(result, "pathway", "pathwayScores")
  expect_equal(dim(scores), c(1L, 6L))
  expect_equal(as.numeric(scores), rep(0, 6))
  expect_identical(rownames(scores), "P1")
  expect_error(createPathwayObject(result, assay = "ramp", slot = "missing",
                                    database = pathwayDb), "not found")
})
